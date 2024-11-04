import base64
import bz2
import gzip
import json
import os
import pathlib
import pickle
from collections import Counter

import numpy as np
import pandas as pd
import pysam

from metaphlan.utils import info, error
from . import utils
from .genotype_resolving import calculate_strain_base_probabilities, resolve_strains
from .get_strainphlan_markers import get_strainphlan_markers
from .linkage import calculate_linkage, linkage_merging
from .model_fitting import fit_model
from .prepare_sam import prepare_sam
from .snp_calling import run_pileup, filter_loci


def try_load_bam(bam_path):
    try:
        return pysam.AlignmentFile(bam_path)
    except ValueError as e:
        if 'file has no sequences defined' in str(e):
            error(f'There are no reads in the BAM file {bam_path}, skipping the sample')
            return None
        else:
            raise e


def step_pileup(sample_path, bam_path, bai_path, output_dir, config):
    info('Preparing the BAM file')
    prepare_sam(sample_path, bam_path, output_dir, config)

    bai_path.unlink(missing_ok=True)
    pysam.index(str(bam_path))

    sam_file = try_load_bam(bam_path)

    if sam_file is None:
        return None, None, None, None

    info('Running pileup')
    base_frequencies, error_rates, marker_to_length = run_pileup(sam_file, config)

    return sam_file, base_frequencies, error_rates, marker_to_length


def save_pileup(pileup_path, base_frequencies, error_rates):
    pileup_path_tmp = pileup_path.with_name(pileup_path.name + '.tmp')
    with gzip.open(pileup_path_tmp, 'wb') as f:
        for m in base_frequencies.keys():
            bfs = {b: base64.b64encode(a) for b, a in base_frequencies[m].items()}
            er = base64.b64encode(error_rates[m])
            line = b'\t'.join([m.encode()] + [bfs[b] for b in 'ACTG'] + [er]) + b'\n'
            f.write(line)

    pileup_path_tmp.rename(pileup_path)


def load_pileup(pileup_path):
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


    return base_frequencies, error_rates, marker_to_length


def step_reconstructed_markers(output_dir_linkages, config, sam_file, base_frequencies, error_rates, marker_to_length,
                               marker_to_clade):
    output_dir_linkages.mkdir(parents=True, exist_ok=True)

    markers_without_clade = [m for m in base_frequencies.keys() if m not in marker_to_clade.keys()]
    if len(markers_without_clade) > 0:
        error(f'There are {len(markers_without_clade)} markers without a clade:'
              f' {utils.shorten_text(",".join(markers_without_clade))}')
        return

    all_markers = sorted(base_frequencies.keys())
    marker_to_clade = {m: marker_to_clade[m] for m in all_markers}
    all_sgbs = sorted(set(marker_to_clade.values()))
    clade_to_markers = pd.Series(marker_to_clade).groupby(marker_to_clade).groups

    result_rows = {}
    strain_resolved_markers = []
    for sgb in all_sgbs:
        loci_rows = []
        for m in clade_to_markers[sgb]:
            for pos in range(marker_to_length[m]):
                bfs = {b: base_frequencies[m][b][pos] for b in 'ACTG'}
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
            continue

        df_loci_sgb = pd.DataFrame(loci_rows)
        marker_to_length_sgb = {m: marker_to_length[m] for m in
                                df_loci_sgb['marker'].unique()}  # downsize the possibly huge dict
        result_row, data, df_loci_sgb, df_loci_sgb_polyallelic = filter_loci(df_loci_sgb, marker_to_length_sgb,
                                                                             config)
        df_loci_biallelic_significant = df_loci_sgb_polyallelic.query('biallelic_significant')
        sgb_linkage = calculate_linkage(df_loci_biallelic_significant, sam_file,
                                        config)  # (m1, p1, m2, p2) => (b1 + b2) => count
        merging_results_before, merging_results = linkage_merging(df_loci_biallelic_significant, sgb_linkage,
                                                                  config)
        sgb_nps = merging_results[0]
        fit_model(df_loci_sgb, result_row)
        df_loci_sgb = calculate_strain_base_probabilities(df_loci_sgb, sgb_nps, result_row, config)
        consensuses_major, consensuses_minor, qualities_major, qualities_minor = resolve_strains(df_loci_sgb,
                                                                                                 marker_to_length_sgb)
        strain_resolved_markers_sgb = [{
            'sgb_id': sgb,
            'marker': m,
            'sequence_maj': consensuses_major[m].decode(),
            'sequence_min': consensuses_minor[m].decode(),
            'qualities_maj': utils.qualities_to_phred(qualities_major[m]),
            'qualities_min': utils.qualities_to_phred(qualities_minor[m]),
        } for m in consensuses_major.keys()]


        with open(output_dir_linkages / f'mrb_{sgb}.pic', 'wb') as f:
            pickle.dump(merging_results_before, f)

        result_rows[sgb] = result_row
        strain_resolved_markers.extend(strain_resolved_markers_sgb)


    # Aggregate output
    df_results_1 = pd.DataFrame.from_dict(result_rows, orient='index')
    per_sgb_info, consensuses_maj, consensuses_min = get_strainphlan_markers(strain_resolved_markers, df_results_1,
                                                                             output_dir_linkages, config)

    df_results_2 = pd.DataFrame.from_dict(per_sgb_info, orient='index')

    df_results = df_results_1.join(df_results_2)

    return df_results, consensuses_maj, consensuses_min


def save_reconstructed_markers(output_results, output_major, output_minor, mp_version, df_results, consensuses_maj,
                               consensuses_min):
    with bz2.open(output_minor, 'wt') as f:
        json.dump({'database_name': mp_version, 'consensus_markers': consensuses_min}, f, indent=2)

    with bz2.open(output_major, 'wt') as f:
        json.dump({'database_name': mp_version, 'consensus_markers': consensuses_maj}, f, indent=2)

    df_results.to_csv(output_results, sep='\t')



def run(sample_path, output_dir, config, target, save_bam_file, reuse, mp_version, marker_to_clade):
    if sample_path.suffix == '.bz2':
        sample_name = sample_path.with_suffix('').stem
    else:
        sample_name = sample_path.stem

    output_dir = output_dir / sample_name

    pileup_path = output_dir / 'pileup.tsv.gz'
    bam_path = output_dir / f'{sample_name}.sorted.bam'
    bai_path = bam_path.with_suffix(bam_path.suffix + '.bai')
    output_dir_linkages = output_dir / 'per_sgb_linkages'
    output_results = output_dir / 'results.tsv'
    output_major = output_dir / f'{sample_name}_major.json.bz2'
    output_minor = output_dir / f'{sample_name}_minor.json.bz2'

    info(f'Starting sample {sample_name}')
    output_dir.mkdir(parents=True, exist_ok=True)

    if save_bam_file:
        files_to_remove = []
    else:
        files_to_remove = [bam_path, bai_path]
    dirs_to_remove = []


    pileup_target_files = [pileup_path]  # files to consider the step to be satisfied
    pileup_output_files = [pileup_path, bam_path, bai_path]  # files to be able to reuse this step
    rm_target_files = [output_results, output_major, output_minor]  # files to consider the step to be satisfied

    if target == 'pileup':
        if reuse in ['all', 'pileup'] and all(map(pathlib.Path.is_file, pileup_target_files)):
            info(f'Pileup for sample {sample_name} is already present, nothing to do')
            return True

        sam_file, base_frequencies, error_rates, marker_to_length = step_pileup(sample_path, bam_path, bai_path,
                                                                                output_dir, config)
        if sam_file is None:
            return False
        save_pileup(pileup_path, base_frequencies, error_rates)

    elif target == 'reconstructed_markers':
        if reuse == 'all' and all(map(pathlib.Path.is_file, rm_target_files)):
            info(f'Reconstructed markers for sample {sample_name} are already present, nothing to do')
            return True

        if reuse in ['all', 'pileup'] and all(map(pathlib.Path.is_file, pileup_output_files)):
            info(f'Loading existing pileup for sample {sample_name}')
            base_frequencies, error_rates, marker_to_length = load_pileup(pileup_path)
            sam_file = try_load_bam(bam_path)
        else:
            info(f'Calculating pileup for sample {sample_name}')
            sam_file, base_frequencies, error_rates, marker_to_length = step_pileup(sample_path, bam_path, bai_path,
                                                                                    output_dir, config)
        if sam_file is None:
            return False

        info(f'Running marker reconstruction for sample {sample_name}')
        df_results, consensuses_maj, consensuses_min = \
            step_reconstructed_markers(output_dir_linkages, config, sam_file, base_frequencies, error_rates,
                                       marker_to_length, marker_to_clade)

        save_reconstructed_markers(output_results, output_major, output_minor, mp_version, df_results, consensuses_maj,
                                   consensuses_min)

        for sgb in df_results.index:
            files_to_remove.append(output_dir_linkages / f'mrb_{sgb}.pic')
        dirs_to_remove.append(output_dir_linkages)


    info(f'Cleaning intermediate files for sample {sample_name}')
    for f in files_to_remove:
        os.unlink(f)

    for d in dirs_to_remove:
        os.rmdir(d)



    info(f'Finished sample {sample_name}')
    return True
