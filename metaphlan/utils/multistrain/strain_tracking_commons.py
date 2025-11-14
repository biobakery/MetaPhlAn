import json
import zipfile
from collections import Counter

import numpy as np

from .utils import MetaphlanDBInfo
from .. import error, warning


def process_samples_argument(samples_list):
    sample_name_to_ac_path = {}
    for ac_path in samples_list:
        suffix = '.allele_counts.zip'
        if not ac_path.name.endswith(suffix):
            error(f'The path {ac_path} does not seem to be allele counts file, skipping it.')
            continue
        sample_name = ac_path.name[:-len(suffix)]
        if not ac_path.exists():
            error(f"The allele counts file does not exist for sample {sample_name}, skipping it.")
            continue
        # sample_dirs.append(sample_dir)
        if sample_name in sample_name_to_ac_path:
            error(f'Duplicated sample name {sample_name}, exiting')
            exit(1)
        sample_name_to_ac_path[sample_name] = ac_path

    if len(sample_name_to_ac_path) == 0:
        error("No samples to run on, exiting")
        return


def load_sample(config, ac_path, target_clade, mp_db_info):
    """

    :param dict config:
    :param pathlib.Path ac_path:
    :param str|None target_clade:
    :param MetaphlanDBInfo mp_db_info:
    :return:
    """

    if not ac_path.exists():
        error(f'Allele counts file does not exist, skipping the sample {ac_path}')
        return None

    clades_all = []
    sample_allele_counts = {}
    sample_marker_breadths = {}
    with zipfile.ZipFile(ac_path, 'r') as f:
        for fi in f.infolist():
            if fi.filename == 'metadata.json':
                with f.open(fi, 'r') as fm:
                    metadata = json.load(fm)
                if 'database_name' in metadata and metadata['database_name'] != mp_db_info.db_name:
                    raise Exception(f'Sample {ac_path.name} has been generated with a different MetaPhlAn database'
                                    f' {metadata["database_name"]} and not the one provided {mp_db_info.db_name}')
            elif fi.filename.endswith('.npy'):
                m = fi.filename[:-len('.npy')]
                clade = mp_db_info.marker_to_clade[m]
                if (target_clade is not None) and (target_clade != clade):
                    continue
                clades_all.append(clade)

                with f.open(fi, 'r') as fm:
                    sample_acs = np.load(fm)
                    sample_allele_counts[m] = sample_acs
                    sample_marker_breadths[m] = (sample_acs.sum(axis=0) > 0).mean()
            else:
                warning(f'Ignoring unrecognized file in the allele counts zip file {fi.filename}')

    clades_to_n_markers = Counter(clades_all)  # clade to number of marker with non-zero coverage

    markers_filtered = sample_allele_counts.keys()
    markers_non_quasi = [m for m in markers_filtered if len(mp_db_info.marker_to_ext[m]) == 0]
    markers_quasi = [m for m in markers_filtered if len(mp_db_info.marker_to_ext[m]) > 0]
    # if any external SGB has >= quasi_marker_frac non-zero markers ==> discard
    markers_quasi_filtered = [m for m in markers_quasi
                              if not any(clades_to_n_markers[ext_sgb] / mp_db_info.clade_to_n_markers[ext_sgb] >=
                                         config['quasi_marker_frac']
                                         for ext_sgb in mp_db_info.marker_to_ext[m])]

    markers_filtered = markers_non_quasi + markers_quasi_filtered
    markers_filtered = [m for m in markers_filtered if sample_marker_breadths[m] >= config['min_breadth']]

    clades_to_n_filtered_markers = Counter([mp_db_info.marker_to_clade[m] for m in markers_filtered])

    clades_present = set([sgb_id for sgb_id, c in clades_to_n_filtered_markers.items()
                          if c >= config['min_markers_abs']
                          and c / mp_db_info.clade_to_n_markers[sgb_id] >= config['min_markers_rel']])
    markers_filtered_2 = set([m for m in markers_filtered if mp_db_info.marker_to_clade[m] in clades_present])

    sample_allele_counts = {m: c for m, c in sample_allele_counts.items() if m in markers_filtered_2}

    return clades_present, sample_allele_counts


def load_sample_parallel(arg):
    config = load_sample_parallel.config
    mp_db_info = load_sample_parallel.mp_db_info

    sample_name, allele_counts_path = arg

    return sample_name, load_sample(config, allele_counts_path, None, mp_db_info)


def load_sample_parallel_only_clades(arg):
    sample_name, allele_counts_path = arg

    res = load_sample_parallel(allele_counts_path)
    if res is None:
        return None

    clades_present, _ = res

    return sample_name, clades_present
