import base64
import bz2
import gzip
import json
import os
import pickle
from collections import Counter

import numpy as np
import pandas as pd
import pysam
from tqdm.auto import tqdm

from . import utils
from metaphlan.utils import info, warning, error
from .get_strainphlan_markers import get_strainphlan_markers
from .genotype_resolving import calculate_strain_base_probabilities, resolve_strains
from .linkage import calculate_linkage, linkage_merging
from .model_fitting import fit_model
from .prepare_sam import prepare_sam
from .snp_calling import run_pileup, filter_loci, estimate_error_rate_pcr_duplicates, \
    estimate_error_rate_overlaps, reuse_error_rate, save_error_rate


def run(sample_path, output_dir, config, save_intermediate_files, reuse, mp_version, marker_to_clade):
    pileup_path = output_dir / 'pileup.tsv.gz'

    if sample_path.suffix == '.bz2':
        sample_name = sample_path.with_suffix('').stem
    else:
        sample_name = sample_path.stem

    info(f'Starting sample {sample_name}')

    if reuse == 'all' and (output_dir / f'{sample_name}_major.json.bz2').exists() and pileup_path.exists():
        info(f'{sample_name} Output already exists, exiting.')
        return

    info(f'{sample_name} started.')

    output_dir.mkdir(parents=True, exist_ok=True)
    bam_path = output_dir / '{}.sorted.bam'.format(sample_name)
    bai_path = bam_path.with_suffix(bam_path.suffix + '.bai')

    if not bam_path.exists() or reuse == 'none' or (reuse != 'none' and not bai_path.exists()):
        info('Preparing the BAM file')
        bam_path.unlink(missing_ok=True)
        bai_path.unlink(missing_ok=True)
        prepare_sam(sample_path, bam_path, output_dir, config)
    else:
        info('Reusing the already prepared BAM file')

    # BAM file index
    if not bai_path.exists() or reuse == 'none':
        info('Indexing the BAM file')
        bai_path.unlink(missing_ok=True)
        pysam.index(str(bam_path))
    else:
        info('The BAM file is already indexed')

    info('Loading the BAM file')
    try:
        sam_file = pysam.AlignmentFile(bam_path)
    except ValueError as e:
        if 'file has no sequences defined' in str(e):
            warning('Empty bam file, exiting')
            return
        else:
            raise e


    if config['error_rate'] == 'estimate':
        error_rate_info_file = output_dir / 'error_rate.info'
        if error_rate_info_file.exists() and reuse == 'all':
            matches_1, matches_2, mismatches_1, mismatches_2 = reuse_error_rate(error_rate_info_file)
        else:
            info('Estimating error rate')
            matches_1, mismatches_1, qual_c = estimate_error_rate_pcr_duplicates(sam_file, config)
            matches_2, mismatches_2 = estimate_error_rate_overlaps(sam_file, config)
            save_error_rate(matches_1, matches_2, mismatches_1, mismatches_2, qual_c, config, error_rate_info_file)

        matches = matches_1 + matches_2
        mismatches = mismatches_1 + mismatches_2
        n_positions = matches + mismatches
        if n_positions > config['error_rate_estimation_positions']:
            err_rate = mismatches / n_positions
            e_phred = -10 * np.log10(err_rate)
            info(f'Using estimated sequencing error rate {err_rate:.1e} (empirical PHRED {e_phred:.1f})')
        else:
            err_rate = config['error_rate_fallback']
            if err_rate == 'phred':
                info('Using fallback error rate calculated per base from base qualities')
            else:
                e_phred = -10 * np.log10(err_rate)
                info(f'Not enough positions to estimated sequencing error rate ({n_positions}). '
                     f'Using fallback error rate {err_rate:.1e} (empirical PHRED {e_phred:.1f})')
    else:
        err_rate = config['error_rate']

    # >> Run pileup <<
    ACTG = 'ACTG'
    pileup_done_path = output_dir / 'pileup.done'
    if pileup_done_path.exists() and reuse in ['pileup', 'all']:
        info('Reusing existing pileup')
        with gzip.open(pileup_path, 'rb') as f:
            base_frequencies = {}
            error_rates = {}
            marker_to_length = {}
            for line in f:
                bf = {}
                m, bf['A'], bf['C'], bf['T'], bf['G'], er = line.strip().split(b'\t')
                m = m.decode()
                base_frequencies[m] = {b: np.frombuffer(base64.decodebytes(x), dtype=np.uint16) for b, x in bf.items()}
                error_rates[m] = np.frombuffer(base64.decodebytes(er), dtype=np.float64)
                marker_to_length[m] = len(error_rates[m])
    else:
        info('Running the pileup')
        base_frequencies, error_rates, marker_to_length = run_pileup(sam_file, config, err_rate)
        info('Saving the pileup')
        with gzip.open(pileup_path, 'wb') as f:
            for m in base_frequencies.keys():
                bfs = {b: base64.b64encode(a) for b, a in base_frequencies[m].items()}
                er = base64.b64encode(error_rates[m])
                line = b'\t'.join([m.encode()] + [bfs[b] for b in ACTG] + [er]) + b'\n'
                f.write(line)
        pileup_done_path.touch()


    markers_without_clade = [m for m in base_frequencies.keys() if m not in marker_to_clade.keys()]
    if len(markers_without_clade) > 0:
        error(f'There are {len(markers_without_clade)} markers without a clade:'
              f' {utils.shorten_text(",".join(markers_without_clade))}')
        return

    all_markers = sorted(base_frequencies.keys())
    marker_to_clade = {m: marker_to_clade[m] for m in all_markers}
    all_sgbs = sorted(set(marker_to_clade.values()))
    clade_to_markers = pd.Series(marker_to_clade).groupby(marker_to_clade).groups

    output_multistrain = output_dir / f'{sample_name}_multistrain.json.gz'
    if reuse == 'all' and output_multistrain.exists():
        with gzip.open(output_multistrain, 'rt') as f:
            x = json.load(f)

        database_name = x['database_name']
        if database_name != mp_version:
            error('Reused multistrain file has different database version')

        strain_resolved_markers = x['strain_resolved_markers']

        df_results = pd.read_csv(output_dir / 'results.tsv', sep='\t', index_col=0)
    else:
        info('Running per SGB')
        result_rows = {}
        strain_resolved_markers = []
        for sgb in tqdm(all_sgbs):
            output_dir_sgb = output_dir / 'per_sgb' / sgb
            progress_file = output_dir_sgb / '_progress.done'

            if progress_file.exists() and reuse == 'all':
                if not (output_dir_sgb / 'result.tsv').exists():
                    continue

                result_row = pd.read_csv(output_dir_sgb / 'result.tsv', sep='\t', index_col=0,
                                         dtype={'snp_rate': float}).squeeze('columns').to_dict()

                with gzip.open(output_dir_sgb / f'multistrain.json.gz', 'rt') as f:
                    json_markers = json.load(f)
                    assert json_markers['database_name'] == mp_version
                    strain_resolved_markers_sgb = json_markers['strain_resolved_markers']
            else:
                loci_rows = []
                for m in clade_to_markers[sgb]:
                    for pos in range(marker_to_length[m]):
                        bfs = {b: base_frequencies[m][b][pos] for b in ACTG}
                        bfs = Counter({k: v for k, v in bfs.items() if v > 0})
                        if bfs.total() == 0:
                            continue
                        loci_rows.append({
                            'marker': m,
                            'pos': pos + 1,
                            'error_rate': error_rates[m][pos],
                            'base_frequencies': bfs,
                        })
                if len(loci_rows) == 0:
                    output_dir_sgb.mkdir(parents=True, exist_ok=True)
                    progress_file.touch()
                    continue

                df_loci_sgb = pd.DataFrame(loci_rows)
                marker_to_length_sgb = {m: marker_to_length[m] for m in df_loci_sgb['marker'].unique()}  # downsize the possibly huge dict
                result_row, data, df_loci_sgb, df_loci_sgb_polyallelic = filter_loci(df_loci_sgb, marker_to_length_sgb, config)
                df_loci_biallelic_significant = df_loci_sgb_polyallelic.query('biallelic_significant')
                sgb_linkage = calculate_linkage(df_loci_biallelic_significant, sam_file, config)  # (m1, p1, m2, p2) => (b1 + b2) => count
                merging_results_before, merging_results = linkage_merging(df_loci_biallelic_significant, sgb_linkage, config)
                sgb_nps = merging_results[0]
                fit_model(df_loci_sgb, result_row)
                df_loci_sgb = calculate_strain_base_probabilities(df_loci_sgb, sgb_nps, result_row, config)
                consensuses_major, consensuses_minor, qualities_major, qualities_minor = resolve_strains(df_loci_sgb, marker_to_length_sgb)
                strain_resolved_markers_sgb = [{
                    'sgb_id': sgb,
                    'marker': m,
                    'sequence_maj': consensuses_major[m].decode(),
                    'sequence_min': consensuses_minor[m].decode(),
                    'qualities_maj': utils.qualities_to_phred(qualities_major[m]),
                    'qualities_min': utils.qualities_to_phred(qualities_minor[m]),
                } for m in consensuses_major.keys()]

                output_dir_sgb.mkdir(parents=True, exist_ok=True)
                df_loci_sgb.to_csv(output_dir_sgb / 'loci.tsv.gz', sep='\t', index=False)
                df_loci_sgb_polyallelic.to_csv(output_dir_sgb / 'loci_polyallelic.tsv.gz', sep='\t', index=False)
                pd.Series(result_row).to_csv(output_dir_sgb / 'result.tsv', sep='\t')

                with open(output_dir_sgb / 'linkage.pic', 'wb') as f:
                    pickle.dump(sgb_linkage, f)
                with open(output_dir_sgb / 'merging_results_before.pic', 'wb') as f:
                    pickle.dump(merging_results_before, f)
                with open(output_dir_sgb / 'merging_results_after.pic', 'wb') as f:
                    pickle.dump(merging_results, f)


                with gzip.open(output_dir_sgb / 'multistrain.json.gz', 'wt') as f:
                    json.dump({
                        'database_name': mp_version,
                        'strain_resolved_markers': strain_resolved_markers_sgb,
                    }, f, indent=2)

                progress_file.touch()

            result_rows[sgb] = result_row
            strain_resolved_markers.extend(strain_resolved_markers_sgb)

        info('Aggregating output')

        df_results = pd.DataFrame.from_dict(result_rows, orient='index')
        df_results.to_csv(output_dir / 'results.tsv', sep='\t')

        with gzip.open(output_multistrain, 'wt') as f:
            json.dump({
                'database_name': mp_version,
                'strain_resolved_markers': strain_resolved_markers,
            }, f, indent=2)


    info('Generating marker files for strainphlan')
    per_sgb_info, consensuses_maj, consensuses_min = get_strainphlan_markers(strain_resolved_markers, df_results, output_dir, config)

    info('Writing markers')
    with bz2.open(output_dir / f'{sample_name}_minor.json.bz2', 'wt') as f:
        json.dump({'database_name': mp_version, 'consensus_markers': consensuses_min}, f, indent=2)

    with bz2.open(output_dir / f'{sample_name}_major.json.bz2', 'wt') as f:
        json.dump({'database_name': mp_version, 'consensus_markers': consensuses_maj}, f, indent=2)

    df_results.join(pd.DataFrame.from_dict(per_sgb_info, orient='index')).to_csv(output_dir / 'results.tsv', sep='\t')

    if not save_intermediate_files:
        info('Cleaning intermediate files')
        files_to_remove = [bam_path, bai_path, pileup_done_path]   # keep the jsons and pileup
        dirs_to_remove = []
        for sgb in all_sgbs:
            output_dir_sgb = output_dir / 'per_sgb' / sgb
            files_to_remove.append(output_dir_sgb / 'linkage.pic')
            files_to_remove.append(output_dir_sgb / 'loci.tsv.gz')
            files_to_remove.append(output_dir_sgb / 'loci_polyallelic.tsv.gz')
            files_to_remove.append(output_dir_sgb / 'merging_results_before.pic')
            files_to_remove.append(output_dir_sgb / 'merging_results_after.pic')
            files_to_remove.append(output_dir_sgb / 'multistrain.json.gz')
            files_to_remove.append(output_dir_sgb / '_progress.done')
            files_to_remove.append(output_dir_sgb / 'result.tsv')
            dirs_to_remove.append(output_dir_sgb)

        dirs_to_remove.append(output_dir / 'per_sgb')

        for f in files_to_remove:
            os.unlink(f)

        for d in dirs_to_remove:
            os.rmdir(d)

    info(f'Finished sample {sample_name}')
