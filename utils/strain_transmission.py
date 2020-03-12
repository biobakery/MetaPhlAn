__author__ = ('Aitor Blanco (aitor.blancomiguez@unitn.it), '
             'Mireia Valles-Colomer (mireia.vallescolomer@unitn.it)')
__version__ = '3.0'
__date__ = '21 Feb 2020'

import os, time, sys
import argparse as ap
from .utils import openw, info, error
from pyphlan import PpaTree, dist_matrix


DISTRIBUTION_THRESHOLD = 0.01

"""
Reads and parses the command line arguments of the script.

:returns: the parsed arguments
"""
def read_params():
    p = ap.ArgumentParser(description="")
    p.add_argument('-t', '--tree', type=str, default=None,
                   help="The input tree file")
    p.add_argument('-d', '--dist', type=str, default=None,
                   help="The input PhyPhlAn pairwise distances file")
    p.add_argument('-m', '--metadata', type=str, default=None,
                   help="The input metadata")
    p.add_argument('-o', '--output_dir', type=str, default=None,
                   help="The output directory")
    p.add_argument('--save_dist', action='store_true',
                   help="[Optional] Save the PhyPhlAn pairwise distances file")
    p.add_argument('--csv', action='store_true',
                   help="[Optional] Export the results as CSV")
    p.add_argument('--vertical', action='store_true',
                   help="[Optional] Detect only vertical transmission events")
    p.add_argument('--horizontal', action='store_true',
                   help="[Optional] Detect only horizontal transmission events")
    p.add_argument('--restrictive', action='store_true',
                   help="[Optional] Use only one timepoint of one individual per family/house in the transmission threshold selection")
    p.add_argument('--permissive', action='store_true',
                   help="[Optional] Use all the data in the transmission threshold selection")
    p.add_argument('--threshold', type=float, default=DISTRIBUTION_THRESHOLD,
                   help="[Optional] A custom distribution threshold value. Default: " + str(DISTRIBUTION_THRESHOLD))
    p.add_argument('--projects', type=str, 
                   nargs='+', default=[],
                   help="[Optional] The specific projects in which detect the transmission")
    
    return p.parse_args()


"""
Checks the mandatory command line arguments of the script.

"""
def check_params(args):
    if not args.tree and not args.dist:
        error('-t (or --tree) or -d (or --distances) must be specified', exit=True, 
            init_new_line=True)
    if not args.output_dir:
        error('-o (or --output_dir) must be specified', exit=True, 
            init_new_line=True)
    if not args.metadata:
        error('-m (or --metadata) must be specified', exit=True, 
            init_new_line=True)
    if args.vertical and args.horizontal:
        error('--vertical and --horizontal cannot be specified at the same time',
            exit=True, init_new_line=True)
    if args.permissive and args.restrictive:
        error('--permissive and --restrictive cannot be specified at the same time',
            exit=True, init_new_line=True)
    if args.tree and args.dist:
        error('-t (or --tree) and -d (or --distances) cannot be specified at the same time',
            exit=True, init_new_line=True)


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

            project = line[1]
            house = line[3]
            family = line[4]
            role = line[5]
            subject = line[2]
            timepoint = line[6]
            sample = line[0]

            if project not in info:
                info[project] = dict()
            if house not in info[project]:
                info[project][house] = dict()
            if family not in info[project][house]:
                info[project][house][family] = {'M':{}, 'I':{}, 'U':{}}
            if subject not in info[project][house][family][role]:
                info[project][house][family][role][subject] = dict()
            if timepoint not in info[project][house][family][role][subject]:
                info[project][house][family][role][subject][timepoint] = sample

            samples[sample] = [project, house, family, role, subject, timepoint]

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
:param transmission_type: the type of transmission events to report
:param training_type: the type of the threshold selection
:returns: the training nodes as a list and the samples to metadata dictionary
"""
def get_training_nodes(nodes, metadata, transmission_type, training_type):
    metadata_info, metadata_samples = get_metadata_info(metadata)
    training_nodes = dict()

    for project in metadata_info:
        for house in metadata_info[project]:
            for family in metadata_info[project][house]:
                if training_type == 'r':
                    role = 'M' if bool(metadata_info[project][house][family]['M']) else 'I' if bool(metadata_info[project][house][family]['I']) else 'U'
                    for subject in metadata_info[project][house][family][role]:
                        for timepoint in metadata_info[project][house][family][role][subject]:
                            sample = metadata_info[project][house][family][role][subject][timepoint]
                            training_nodes[sample] = metadata_samples[sample]
                            break
                        break
                    if transmission_type in ['a', 'h']:
                        break
                else:
                    for role in metadata_info[project][house][family]:
                        for subject in metadata_info[project][house][family][role]:
                            for timepoint in metadata_info[project][house][family][role][subject]:
                                sample = metadata_info[project][house][family][role][subject][timepoint]
                                training_nodes[sample] = metadata_samples[sample]
                                if training_type == 'p':
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
def get_training_distances(training_nodes, pairwise_distances, transmission_type):
    training_distances = list()
    for pair in pairwise_distances:
        if pair["1"] in training_nodes and pair["2"] in training_nodes:
            if not training_nodes[pair['1']] or not training_nodes[pair['2']]:
                training_distances.append(pair)
            elif transmission_type == 'v' and not training_nodes[pair['1']][2] == training_nodes[pair['2']][2]:
                training_distances.append(pair)
            elif transmission_type == 'h' and not training_nodes[pair['1']][1] == training_nodes[pair['2']][1]:
                training_distances.append(pair)
            elif (not training_nodes[pair['1']][1] == training_nodes[pair['2']][1]) and (not training_nodes[pair['1']][2] == training_nodes[pair['2']][2]):
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
:param transmission_projects: the projects in which detect the transmission
:param metadata_samples: the sample to metadata dictionary
:param threshold: the calculated threshold
:param transmission_type: the type of transmission events to report
:returns: the transmission events as a list
"""
def get_transmission_events(pairwise_distances, transmission_projects, metadata_samples, threshold, transmission_type):
    transmission_events = list()
    for pair in pairwise_distances:
        if pair['1'] in metadata_samples and pair['2'] in metadata_samples:
            if not(transmission_projects) or (metadata_samples[pair['1']][0] in transmission_projects and metadata_samples[pair['2']][0] in transmission_projects):
                if not metadata_samples[pair['1']][4] == metadata_samples[pair['2']][4] and float(pair['dist']) <= threshold:
                    if transmission_type == 'v' and metadata_samples[pair['1']][2] == metadata_samples[pair['2']][2]:
                        if sorted([metadata_samples[pair['1']][3], metadata_samples[pair['2']][3]]) == ['I', 'M']:
                            transmission_events.append(pair)
                    elif transmission_type == 'h' and metadata_samples[pair['1']][1] == metadata_samples[pair['2']][1]:   
                        if not (sorted([metadata_samples[pair['1']][3], metadata_samples[pair['2']][3]]) == ['I', 'M'] and metadata_samples[pair['1']][2] == metadata_samples[pair['2']][2]):
                            transmission_events.append(pair)
                    elif transmission_type == 'a' and (metadata_samples[pair['1']][1] == metadata_samples[pair['2']][1] or metadata_samples[pair['1']][2] == metadata_samples[pair['2']][2]):   
                        transmission_events.append(pair)

    return transmission_events


"""
Writes the detected transmission events to CSV

:param transmission_events: the detected transmission events
:param output_dir: the output directory to store the results
:param metadata_samples: the sample to metadata dictionary
"""
def write_transmission_events_as_csv(transmission_events, output_dir, metadata_samples):
    with open(os.path.join(output_dir, "transmission_events.csv"), 'w') as report:
        report.write('project\tsample_1\tsample_2\ttransmission_type\tpairwise_distance\tindividual_1\tindividual_2\ttimepoint_1\ttimepoint_2\trole_1\trole_2\tfamily_1\tfamily_2\thouse_1\thouse_2\n')
        for event in transmission_events:
            if sorted([metadata_samples[event['1']][3], metadata_samples[event['2']][3]]) == ['I', 'M'] and metadata_samples[event['1']][2] == metadata_samples[event['2']][2]:  
                t = 'vertical'
                if metadata_samples[event['1']][3] == 'M':
                    e = [event['1'], event['2']]
                else:
                    e = [event['2'], event['1']]
            else:
                t = 'horizontal'
                e = [event['1'], event['2']]
            report.write(metadata_samples[e[0]][0]+ "\t" + e[0] + "\t" + e[1] + '\t' + t + '\t' + event['dist'] + '\t' + metadata_samples[e[0]][4] + '\t' + metadata_samples[e[1]][4] + '\t' + metadata_samples[e[0]][5] + '\t'  +metadata_samples[e[1]][5] + '\t' + metadata_samples[e[0]][3] + '\t'  + metadata_samples[e[1]][3] + '\t' + metadata_samples[e[0]][2] + '\t'  +metadata_samples[e[1]][2] + '\t' + metadata_samples[e[0]][1] + '\t'  + metadata_samples[e[1]][1] + "\n") 


"""
Writes the detected transmission events to file

:param transmission_events: the detected transmission events
:param threshold: the calculated threshold
:param output_dir: the output directory to store the results
:param metadata_samples: the sample to metadata dictionary
:param output_csv: whether export the results as CSV
"""
def write_transmission_events(transmission_events, threshold, output_dir, metadata_samples, output_csv):
    with open(os.path.join(output_dir, "transmission_events.info"), 'w') as report:
        report.write("Selected strain-transmission threshold: "+str(threshold)+"\n")
        report.write("Number of transmission events: "+str(len(transmission_events))+"\n\n")
        for event in transmission_events:
            if sorted([metadata_samples[event['1']][3], metadata_samples[event['2']][3]]) == ['I', 'M'] and metadata_samples[event['1']][2] == metadata_samples[event['2']][2]:
                if metadata_samples[event['1']][3] == 'M':
                    report.write(event['1']+" -> "+event['2']+"\n")
                else:
                    report.write(event['2']+" -> "+event['1']+"\n")
            else:
               report.write(event['1']+" <-> "+event['2']+"\n") 
    if output_csv:
        write_transmission_events_as_csv(transmission_events, output_dir, metadata_samples)


"""
Identifies transmission events in phylogenetic trees

:param tree: the input Newick tree
:param metadata: the metadata file
:param transmission_projects: the projects in which detect the transmission
:param distr_threshold: the distribution threshold
:param transmission_type: the type of transmission events to report
:param save_dist: whether to save the pairwise distances file
:param training_type: the type of the threshold selection
:param output_csv: whether export the results as CSV
:param output_dir: the output directory to store the results
"""
def strain_transmission(tree, distances_file, metadata, transmission_projects, distr_threshold, transmission_type, save_dist, training_type, output_csv, output_dir):
    normalise = True
    matrix = False
    if not distances_file:
        distances_file = tree+".dist"
        tree_pairwisedist(tree, normalise, matrix, os.path.join(output_dir, distances_file))
    pairwise_distances = parse_distances(os.path.join(output_dir, distances_file)) 
    if not distances_file and not save_dist:
        os.remove(os.path.join(output_dir, distances_file))

    nodes = get_nodes(pairwise_distances)
    
    training_nodes, metadata_samples = get_training_nodes(nodes, metadata, transmission_type, training_type)    
    training_distances = get_training_distances(training_nodes, pairwise_distances, transmission_type)

    threshold = get_threshold(training_distances, distr_threshold)

    transmission_events = get_transmission_events(pairwise_distances, transmission_projects, metadata_samples, threshold, transmission_type)
    write_transmission_events(transmission_events, threshold, output_dir, metadata_samples, output_csv)


"""
Main call

:param tree: the input Newick tree
:param metadata: the metadata file
:param projects: the projects in which detect the transmission
:param distr_threshold: the distribution threshold
:param horizontal: whether to detect only horizontal transmission events
:param vertical: whether to detect only vertical transmission events
:param restrictive: whether to be restrictive when threshold selection
:param permissive: whether to be permissive when threshold selection
:param save_dist: whether to save the pairwise distances file
:param csv: whether export the results as CSV
:param output_dir: the output directory to store the results
"""
if __name__ == "__main__":
    t0 = time.time()
    args = read_params()
    info("Start execution")
    check_params(args)

    transmission = 'v' if args.vertical else 'h' if args.horizontal else 'a'
    training = 'r' if args.restrictive else 'p' if args.permissive else 'n'

    strain_transmission(args.tree, args.dist, args.metadata, args.projects, args.threshold, transmission, args.save_dist, training, args.csv, args.output_dir)

    exec_time = time.time() - t0
    info("Finish execution ("+str(round(exec_time, 2))+" seconds)\n", 
        init_new_line=True)
