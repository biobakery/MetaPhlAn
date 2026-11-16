import itertools as it
from collections import defaultdict, Counter

import networkx as nx
import numpy as np

from ...utils import error

ACTG = 'ACTG'

# ==============================================================================================
# Read linkage between significant biallelic sites
# ==============================================================================================
def calculate_linkage(df_loci_sgb, sam_file, config):
    """(marker1, pos1, marker2, pos2) -> {base1 + base2: read names}, from reads (and mate pairs) covering both
    sites with base quality >= min_base_quality. Only the significant sites inside each read are inspected."""
    sgb_linkage = defaultdict(lambda: defaultdict(set))
    sgb_markers = df_loci_sgb['marker'].unique()
    min_bq = config['min_base_quality']

    # significant positions per marker (1-based, sorted)
    site_pos = {m: np.sort(df_loci_sgb.loc[df_loci_sgb['marker'] == m, 'pos'].to_numpy(dtype=np.int64))
                for m in sgb_markers}

    paired_reads = defaultdict(list)
    for marker in sgb_markers:
        for read in sam_file.fetch(marker):
            if read.is_unmapped | read.is_qcfail | read.is_secondary:
                continue
            if config['filter_duplicates'] and read.is_duplicate:
                continue
            if read.mapping_quality < config['min_mapping_quality']:
                continue
            paired_reads[read.qname].append(read)

    for read_pair_name, reads in paired_reads.items():
        if len(reads) > 2:
            error(f"Something wrong with paired-end, more than two ({len(reads)}) mates found with name {read_pair_name}")
            continue
        if len(reads) == 2 and (reads[0].reference_name == reads[1].reference_name) \
                and not (reads[0].is_forward ^ reads[1].is_forward):
            continue

        # per mate: {(marker, pos): base} at significant positions with base quality >= threshold
        mate_obs = []
        for read in reads:
            marker = read.reference_name
            sp_ = site_pos[marker]
            ct = read.cigartuples
            qual = read.query_qualities
            obs = {}
            if ct is not None and len(ct) == 1 and ct[0][0] == 0 and qual is not None \
                    and read.query_length == ct[0][1]:
                # one M operation with SEQ and QUAL
                start = read.reference_start
                L = ct[0][1]
                lo = np.searchsorted(sp_, start + 1, 'left')
                hi = np.searchsorted(sp_, start + L, 'right')
                if hi > lo:
                    seq = read.query_alignment_sequence
                    for p in sp_[lo:hi]:
                        p = int(p)
                        o = p - 1 - start
                        if qual[o] >= min_bq:
                            obs[(marker, p)] = seq[o]
            else:  # any other read: aligned pairs from pysam, restricted to significant sites
                ref_pos = read.get_reference_positions(full_length=True)
                read_seq = read.query_alignment_sequence
                read_quals = read.query_qualities

                assert len(ref_pos) == len(read_seq)
                assert len(read_quals) == len(read_seq)

                sset = set(int(x) for x in sp_)
                for p, b, q in zip(ref_pos, read_seq, read_quals):
                    if p is not None and q >= min_bq and (p + 1) in sset:
                        obs[(marker, p + 1)] = b
            mate_obs.append(obs)

        covered = {}
        if len(mate_obs) == 1:
            covered = mate_obs[0]
        else:
            a, b = mate_obs
            covered = dict(a)
            for k, v in b.items():
                if k in covered:
                    if covered[k] != v:
                        covered[k] = None  # mates disagree -> drop the site
                else:
                    covered[k] = v
            covered = {k: v for k, v in covered.items() if v is not None}

        if len(covered) < 2:
            continue
        covered_loci = sorted(covered)
        for (m1, p1), (m2, p2) in it.combinations(covered_loci, 2):
            sgb_linkage[(m1, p1, m2, p2)][covered[(m1, p1)] + covered[(m2, p2)]].add(read_pair_name)

    return {k1: {k2: v2 for k2, v2 in v1.items()} for k1, v1 in sgb_linkage.items()}


class Node:
    def __init__(self, marker_pos, genotype):
        self.marker_pos = marker_pos  # array of (marker, pos) pairs
        self.genotype = genotype

    def __repr__(self):
        return '|'.join(m + '-' + str(p) for m, p in self.marker_pos) + self.genotype

    def __str__(self):
        return self.__repr__()


class NodePair:
    def __init__(self, marker_pos, genotype_maj, genotype_min):
        self.marker_pos = marker_pos  # array of (marker, pos) pairs
        self.genotype_maj = genotype_maj
        self.genotype_min = genotype_min

    def get_node_major(self):
        return Node(self.marker_pos, self.genotype_maj)

    def get_node_minor(self):
        return Node(self.marker_pos, self.genotype_min)


def eval_pair(g, np1, np2):
    """

    :param g:
    :param NodePair np1:
    :param NodePair np2:
    :return:
    """
    n1_a = np1.get_node_minor()
    n1_b = np1.get_node_major()
    n2_a = np2.get_node_minor()
    n2_b = np2.get_node_major()
    weights = []
    for n1 in (n1_a, n1_b):
        for n2 in (n2_a, n2_b):
            n1_name = str(n1)
            n2_name = str(n2)
            if n2_name in g[n1_name]:
                w = len(g[n1_name][n2_name]['reads'])
            else:
                w = 0
            weights.append(w)
    w1 = (weights[0], weights[3])
    w2 = (weights[1], weights[2])


    if sum(w1) < sum(w2):  # make sure w1 is along, w2 is across; a is linked with a, b is linked with b
        w1, w2 = w2, w1
        np2 = NodePair(np2.marker_pos, np2.genotype_min, np2.genotype_maj)  # swap min/maj

    return {'np1': np1, 'np2': np2, 'w_along': w1, 'w_across': w2}


def merge(r, g):
    old_np_1 = r['np1']
    old_np_2 = r['np2']
    old_node_name_1a = str(old_np_1.get_node_major())
    old_node_name_1b = str(old_np_1.get_node_minor())
    old_node_name_2a = str(old_np_2.get_node_major())
    old_node_name_2b = str(old_np_2.get_node_minor())
    new_marker_pos = old_np_1.marker_pos + old_np_2.marker_pos
    new_gt_a = old_np_1.genotype_maj + old_np_2.genotype_maj
    new_gt_b = old_np_1.genotype_min + old_np_2.genotype_min
    new_node_name_a = str(Node(new_marker_pos, new_gt_a))
    new_node_name_b = str(Node(new_marker_pos, new_gt_b))
    g.add_node(new_node_name_a)
    g.add_node(new_node_name_b)

    for new_node_name_x, old_node_name_1x, old_node_name_2x in ((new_node_name_a, old_node_name_1a, old_node_name_2a), (new_node_name_b, old_node_name_1b, old_node_name_2b)):
        for other_node in set.union(set(g[old_node_name_1x]), set(g[old_node_name_2x])):
            if other_node in (new_node_name_a, new_node_name_b, old_node_name_1a, old_node_name_1b, old_node_name_2a, old_node_name_2b):
                continue
            reads_x = set()
            if other_node in g[old_node_name_1x]:
                reads_x |= g[old_node_name_1x][other_node]['reads']
            if other_node in g[old_node_name_2x]:
                reads_x |= g[old_node_name_2x][other_node]['reads']
            if len(reads_x) > 0:
                g.add_edge(new_node_name_x, other_node, reads=reads_x)

    g.remove_nodes_from([old_node_name_1a, old_node_name_1b, old_node_name_2a, old_node_name_2b])

    np_merged = NodePair(new_marker_pos, new_gt_a, new_gt_b)
    return np_merged


def generate_node_pairs(df_loci_sgb):
    node_pairs = []   # triplet (list of positions, string of major bases, string of minor bases)
    nodes = {}
    for i, ((marker, pos), row) in enumerate(df_loci_sgb.iterrows()):
        base_frequencies = row['base_frequencies']
        base_frequencies = Counter(dict(zip(ACTG, base_frequencies)))
        b_maj, b_min = [x[0] for x in base_frequencies.most_common(2)]
        node_pair = NodePair([(marker, pos)], b_maj, b_min)
        node_pairs.append(node_pair)
        for j, n in enumerate((node_pair.get_node_minor(), node_pair.get_node_major())):
            node_name = str(n)
            ratio = base_frequencies[n.genotype] / row['base_coverage']
            data = {
                'freq': base_frequencies[n.genotype],
                'base_coverage': row['base_coverage'],
                'ratio': ratio
            }
            nodes[node_name] = data

    return nodes, node_pairs


def generate_graph(nodes, edges):
    g = nx.Graph()
    g.add_nodes_from(nodes)
    for n in nodes:
        for k, v in nodes[n].items():
            g.nodes[n][k] = v  # add attributes to the nodes

    edges_ = []
    for (m1, p1, m2, p2), gt_c in edges.items():
        for (b1, b2), c in gt_c.items():
            n1 = str(Node([(m1, p1)], b1))
            n2 = str(Node([(m2, p2)], b2))
            edges_.append((n1, n2, {'reads': c}))

    g.add_edges_from(edges_)
    return g


def linkage_merging(df_loci_sgb, sgb_linkage, config):
    nodes, node_pairs = generate_node_pairs(df_loci_sgb.set_index(['marker', 'pos']))
    g = generate_graph(nodes, sgb_linkage)

    if not config['linkage_graph_merging']:
        return [node_pairs, g], [node_pairs, g]

    results_before_merging = [node_pairs, g.copy()]
    node_pairs_merged = []

    np1 = None
    for np2 in node_pairs:
        if np1 is None:
            np1 = np2
            continue

        r = eval_pair(g, np1, np2)
        if np1.marker_pos[0][0] == np2.marker_pos[0][0] and sum(r['w_across']) <= config['max_w_across'] \
                and sum(r['w_along']) >= config['min_w_along']:  # only same marker
            np1 = merge(r, g)
        else:
            node_pairs_merged.append(np1)
            np1 = np2
    if np1 is not None:
        node_pairs_merged.append(np1)

    unmerged_to_merged = {}
    for node_pair in node_pairs_merged:
        for m, p in node_pair.marker_pos:
            assert (m, p) not in unmerged_to_merged
            unmerged_to_merged[(m, p)] = node_pair

    cross_pairs = set()
    for k in sgb_linkage.keys():
        m1, p1, m2, p2 = k
        if m1 != m2:
            cross_pairs.add((unmerged_to_merged[(m1, p1)], unmerged_to_merged[(m2, p2)]))

    while len(cross_pairs) > 0:
        np1, np2 = cross_pairs.pop()
        if (np2, np1) in cross_pairs:  # solve the problem of ordering
            cross_pairs.remove((np2, np1))

        r = eval_pair(g, np1, np2)
        if sum(r['w_across']) <= config['max_w_across'] and sum(r['w_along']) >= config['min_w_along']:
            np_merged = merge(r, g)

            for m, p in np_merged.marker_pos:
                unmerged_to_merged[(m, p)] = np_merged

            to_remove = set()
            to_add = set()
            for np1_, np2_ in cross_pairs:
                if np1_ == np1 or np1_ == np2:
                    to_remove.add((np1_, np2_))
                    to_add.add((np_merged, np2_))
                if np2_ == np1 or np2_ == np2:
                    to_remove.add((np1_, np2_))
                    to_add.add((np1_, np_merged))

            if not set(to_remove).issubset(cross_pairs):
                print(np1, np2)
                print(cross_pairs)
                print(to_remove)
                print(to_add)
                assert False

            cross_pairs = cross_pairs.difference(to_remove).union(to_add)

    node_pairs_after_cross_merging = list(set(unmerged_to_merged.values()))

    results = [node_pairs_after_cross_merging, g]

    return results_before_merging, results
