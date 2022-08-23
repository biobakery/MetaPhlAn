__author__ = ('Aitor Blanco (aitor.blancomiguez@unitn.it), '
             'Mireia Valles-Colomer (mireia.vallescolomer@unitn.it)')
__version__ = '4.0.0'
__date__ = '22 Aug 2022'

import os, time, sys
import argparse as ap

try:
    from .util_fun import openw, info, error
    from .pyphlan import PpaTree, dist_matrix
except ImportError:
    from util_fun import openw, info, error
    from pyphlan import PpaTree, dist_matrix


DISTRIBUTION_THRESHOLD = 0.03
metaphlan_script_install_folder = os.path.dirname(os.path.abspath(__file__))
DEFAULT_UTILS_FOLDER = os.path.join(metaphlan_script_install_folder)
DEFAULT_UTILS_FOLDER = os.environ.get('METAPHLAN_DB_DIR', DEFAULT_UTILS_FOLDER)
PRECOMPUTED_FILE = os.path.join(DEFAULT_UTILS_FOLDER, 'VallesColomerM_2022_Nov19_thresholds.tsv')


"""
Reads and parses the command line arguments of the script.

:returns: the parsed arguments
"""
def read_params():
    p = ap.ArgumentParser(description="", formatter_class=ap.ArgumentDefaultsHelpFormatter)
    p.add_argument('-t', '--tree', type=str, default=None,
                   help="The input tree file")
    p.add_argument('-m', '--metadata', type=str, default=None,
                   help="The input metadata")
    p.add_argument('-o', '--output_dir', type=str, default=None,
                   help="The output directory")
    p.add_argument('--sgb_id', type=str, default=None,
                   help="[Optional] If specified, it will use the precomputed transmisison threshold for the specific SGB from the VallesColomerM_2022 study")
    p.add_argument('--threshold', type=float, default=DISTRIBUTION_THRESHOLD,
                   help="[Optional] A custom distribution threshold value")
    p.add_argument('--precomputed_thresholds_file', type=str, default=PRECOMPUTED_FILE,
                   help="[Optional] The file containing the pre-computed thresholds")
    p.add_argument('--save_dist', action='store_true',
                   help="[Optional] Save the PhyPhlAn pairwise distances file")
    
    return p.parse_args()


"""
Checks the mandatory command line arguments of the script.

"""
def check_params(args):
    if not args.tree and not args.dist:
        error('-t (or --tree) must be specified', exit=True, 
            init_new_line=True)
    if not args.output_dir:
        error('-o (or --output_dir) must be specified', exit=True, 
            init_new_line=True)
    elif not os.path.exists(args.output_dir):
        error('The directory {} does not exist'.format(args.output_dir), exit=True, 
            init_new_line=True)
    if not args.metadata:
        error('-m (or --metadata) must be specified', exit=True, 
            init_new_line=True)


"""
Calls the main function on the tree_pairwisedist.py Script of the Pyphlan tool

:param tree: the input Newick tree
:normalise: whether normalise the distances with respect to the total branch length
:matrix: whether report the results as a matrix
:output: the output distances file
"""
def tree_pairwisedist(tree, normalise, matrix, output):
    ppatree = PpaTree(tree)

    dists = dist_matrix(ppatree.tree) 
    tbl = ppatree.tree.total_branch_length() if normalise else 1.0

    with openw(output) as out:
        if matrix:
            keys = sorted(dists.keys())
            out.write( "\t".join(["ID"]+keys) +"\n" )
            for k1 in keys:
                out.write( "\t".join([k1]+[str(dists[k1][k2]/tbl) for k2 in keys]) +"\n" )
        else:
            for k1,v1 in dists.items():
                for k2,v2 in v1.items():
                    if k1 < k2:
                        out.write( "\t".join([k1,k2,str(v2/tbl)]) +"\n" )


"""
Gets all the information about the families from the metadata

:param metadata: the metadata
:returns: the tree about the samples in the metadata and
    the samples to metadata dictionary
"""
def get_metadata_info(metadata):
    info = dict()
    samples = dict()
    with open(metadata, 'r') as metadata_file:
        metadata_file.readline()
        for line in metadata_file:
            line = line.strip().split("\t")

            relation = line[2]
            subject = line[1]
            timepoint = line[3]
            sample = line[0]

            
            if relation not in info:
                info[relation] = dict()
            if subject not in info[relation]:
                info[relation][subject] = dict()
            if timepoint not in info[relation][subject]:
                info[relation][subject][timepoint] = sample

            samples[sample] = [relation, subject, timepoint]

    return info, samples


""""
Parses the pairwise distances from the generated distances file

:param distances_path: the pairwise distances file
:returns: the pairwise distances as a list
"""
def parse_distances(distances_path):
    distances = list()
    with open(distances_path, 'r') as distances_file:
        for line in distances_file:
            line = line.strip().split("\t")
            distances.append({"1": line[0], "2": line[1], "dist": line[2]})
    return distances


"""
Gets the nodes of the tree

:param distances: the pairwise distances from the tree
:returns: the nodes of the tree as a set
"""
def get_nodes(distances):
    nodes = set()
    for line in distances:
        nodes.add(line["1"])
        nodes.add(line["2"])
    return nodes


"""
Gets the training nodes from the metadata

:param nodes: the tree nodes
:param metadata: the metadata file
:returns: the training nodes as a list and the samples to metadata dictionary
"""
def get_training_nodes(nodes, metadata):
    metadata_info, metadata_samples = get_metadata_info(metadata)
    training_nodes = dict()

    for relation in metadata_info:
        for subject in metadata_info[relation]:
            for timepoint in metadata_info[relation][subject]:
                sample = metadata_info[relation][subject][timepoint]
                training_nodes[sample] = metadata_samples[sample]
                break

    for node in nodes:
        if node not in metadata_samples:
            training_nodes[node]=list()

    return training_nodes, metadata_samples


"""
Gets the distances between the training nodes

:param training_nodes: the training nodes
:param pairwise_distances: a list with the pairwise distances
:param transmission_type: the type of transmission events to report
:returns: the training distances as a list
"""
def get_training_distances(training_nodes, pairwise_distances):
    training_distances = list()
    for pair in pairwise_distances:
        if pair["1"] in training_nodes and pair["2"] in training_nodes:
            if not training_nodes[pair['1']] or not training_nodes[pair['2']]:
                training_distances.append(pair)
            elif not training_nodes[pair['1']][0] == training_nodes[pair['2']][0]:
                training_distances.append(pair)

    return training_distances


"""
Gets the threshold defining strain transmission

:param training_distances: the distances between the training nodes
:param distr_threshold: the distribution threshold
:returns: the strain transmission threshold
"""
def get_threshold(training_distances, distr_threshold):
    distances = list()
    for distance in training_distances:
        distances.append(float(distance['dist']))
    distances.sort()
    return distances[int(len(distances)*distr_threshold)]

    
"""
Gets the transmission events using a calculated threshold

:param pairwise_distances: the pairwise distances
:param metadata_samples: the sample to metadata dictionary
:param threshold: the calculated threshold
:returns: the transmission events as a list
"""
def get_transmission_events(pairwise_distances, metadata_samples, threshold):
    transmission_events = list()
    for pair in pairwise_distances:
        if pair['1'] in metadata_samples and pair['2'] in metadata_samples:
            if not metadata_samples[pair['1']][1] == metadata_samples[pair['2']][1] and float(pair['dist']) <= threshold:
                if metadata_samples[pair['1']][0] == metadata_samples[pair['2']][0]:   
                    transmission_events.append(pair)

    return transmission_events


"""
Writes the detected transmission events to file

:param transmission_events: the detected transmission events
:param threshold: the calculated threshold
:param output_dir: the output directory to store the results
:param metadata_samples: the sample to metadata dictionary
"""
def write_transmission_events(transmission_events, threshold, output_dir, metadata_samples):
    with open(os.path.join(output_dir, "transmission_events.info"), 'w') as report:
        report.write("Selected strain-transmission threshold: "+str(threshold)+"\n")
        report.write("Number of transmission events: "+str(len(transmission_events))+"\n\n")
        for event in transmission_events:
            report.write(event['1']+" <-> "+event['2']+"\n") 

"""
Gets the precomputed threshold from VallesColomerM_2022 study

:param sgb_id: the SGB id
:returns: the transmission threshold
"""
def get_precomputed_threshold(sgb_id, precomputed_thresholds_file):
    sgb2thres = dict()
    with open(precomputed_thresholds_file, 'r') as rf:
        rf.readline()
        for line in rf:
            line = line.strip().split('\t')
            sgb2thres[line[0]] = line[1]
    if sgb_id not in sgb2thres:
        error('The SGB specified "{}" has not been precomputed'.format(sgb_id), exit=True, init_new_line=True)
    else:
        return float(sgb2thres[sgb_id])

"""
Identifies transmission events in phylogenetic trees

:param tree: the input Newick tree
:param metadata: the metadata file
:param distr_threshold: the distribution threshold
:param sgb_id: the SGB id
:param precomputed_thresholds_file: the file containing the precomputed thresholds
:param save_dist: whether to save the pairwise distances file
:param output_dir: the output directory to store the results
"""
def strain_transmission(tree, metadata, distr_threshold, sgb_id, precomputed_thresholds_file, save_dist, output_dir):
    normalise = True
    matrix = False
    distances_file = tree+".dist"
    tree_pairwisedist(tree, normalise, matrix, os.path.join(output_dir, distances_file))
    pairwise_distances = parse_distances(os.path.join(output_dir, distances_file)) 
    if not save_dist:
        os.remove(os.path.join(output_dir, distances_file))

    nodes = get_nodes(pairwise_distances)
    training_nodes, metadata_samples = get_training_nodes(nodes, metadata)    
    training_distances = get_training_distances(training_nodes, pairwise_distances)
    if sgb_id is None:
        threshold = get_threshold(training_distances, distr_threshold)
    else:
        threshold = get_threshold(training_distances, get_precomputed_threshold(sgb_id, precomputed_thresholds_file))

    transmission_events = get_transmission_events(pairwise_distances, metadata_samples, threshold)
    write_transmission_events(transmission_events, threshold, output_dir, metadata_samples)


"""
Main call

:param tree: the input Newick tree
:param metadata: the metadata file
:param threshold: the distribution threshold
:param sgb_id: the SGB id
:param precomputed_thresholds_file: the file containing the precomputed thresholds
:param save_dist: whether to save the pairwise distances file
:param output_dir: the output directory to store the results
"""
def main():
    t0 = time.time()
    args = read_params()
    info("Start execution")
    check_params(args)

    strain_transmission(args.tree, args.metadata, args.threshold, args.sgb_id, args.precomputed_thresholds_file, args.save_dist, args.output_dir)

    exec_time = time.time() - t0
    info("Finish execution ("+str(round(exec_time, 2))+" seconds)\n", 
        init_new_line=True)
	
if __name__ == "__main__":
	main()

