#! /usr/bin/env python

from .scripts import PROGRAM_VERSION, PROGRAM_NAME
from .scripts import __file__ as treeshrink_file
from .scripts.sequence_lib import sample_from_list
from .scripts.optimal_filter_lib import TreeFilter
from .scripts.tree_lib import prune_tree, get_taxa,tree_as_newick
from sys import argv, stdout,setrecursionlimit
from math import sqrt,log,exp
from subprocess import check_output,call
import argparse
from dendropy import Tree, TreeList
from os.path import basename, dirname, splitext,realpath,join,normpath,isdir,isfile,exists
from os import mkdir,getcwd,rmdir,listdir
from copy import deepcopy
from shutil import rmtree, copyfile
from .scripts.alignment import CompactAlignment
from .scripts import set_tmp_dir, get_tmp_dir, get_tmp_file
from .scripts.util_lib import minVar_bisect
import re
import random    

def make_dir(dirName):
    if exists(dirName) and isdir(dirName):
        return False
    mkdir(dirName)
    return True

def test_Rlib(libdir):
    filename = get_tmp_file("test_Rlib.txt")    
    with open(filename,'w') as fout:
        for i in range(300):
            fout.write(str(1+random.lognormvariate(0,1)) + "\n")
    try:
        check_output(["Rscript",normpath(join(libdir,"R_scripts","find_threshold_lkernel.R")),libdir,filename,"0.05"]).lstrip().rstrip()[5:]
        return True
    except:
        return False

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i","--indir",required=False,help="The parent input directory where the trees (and alignments) can be found")
    parser.add_argument("-t","--tree",required=False,default="input.tree",help="The name of the input tree/trees. If the input directory is specified (see -i option), each subdirectory under it must contain a tree with this name. Otherwise, all the trees can be included in this one file. Default: input.tree")
    parser.add_argument("-g","--g2sp",required=False,help="The gene-name-to-species-name mapping file")
    parser.add_argument("-a","--alignment",required=False,help="The name of the input alignment; can only be used when the input directory is specified (see -i option). Each subdirectory under it must contain an alignment with this name. Default: input.fasta")
    parser.add_argument("-c","--centroid",required=False,action='store_true',help="Do centroid reroot in preprocessing. Highly recommended for large trees. Default: NO")
    parser.add_argument("-k","--k",required=False,help="The maximum number of leaves that can be removed. Default: auto-select based on the data; see also -s")
    parser.add_argument("-s","--kscaling",required=False,help="If -k not given, we use k=min(n/a,b*sqrt(n)) by default; using this option, you can set the a,b constants; Default: '5,2'")
    parser.add_argument("-q","--quantiles",required=False,help="The quantile(s) to set threshold. Default is 0.05")
    parser.add_argument("-b","--minImpact",required=False,help="Do not remove species on the per-species test if their impact on diameter is less than x%% where x is the given value. Default: 5")
    parser.add_argument("-m","--mode",required=False,help="Filtering mode: 'per-species', 'per-gene', 'all-genes','auto'. Default: auto")
    parser.add_argument("-o","--outdir",required=False,help="Output directory. Default: If the input directory is specified, outputs will be placed in that input directory. Otherwise, a directory with the suffix 'treeshrink' will be created in the same place as the input trees")
    parser.add_argument("-O","--outprefix",default="output",required=False,help="Output name prefix. If the output directory contains some files with the specified prefix, automatically adjusts the prefix (e.g. output --> output1) to avoid overriding. Use --force to force overriding. Default: 'output'")
    parser.add_argument("-f","--force",required=False,action='store_true',help="Force overriding of existing output files.")
    parser.add_argument("-p","--tempdir",required=False,help="Directory to keep temporary files. If specified, the temp files will be kept")
    parser.add_argument("-v","--version",required=False,action='store_true',help="Show TreeShrink version.")
    parser.add_argument("-x","--exceptions",required=False,help="A list of special species that will not be removed in any of the input trees.")

    if len(argv) == 1:
        parser.print_help()
        exit(0)

    args = vars(parser.parse_args())
    
    if args["version"]:
        print(treeshrink.PROGRAM_VERSION)
        exit(0)

    setrecursionlimit(5000)

    print("Launching " + PROGRAM_NAME + " version " + PROGRAM_VERSION)
    print(PROGRAM_NAME + " was called as follow")
    print(" ".join(argv))


    MIN_OCC = 20
    MIN_TREE_NUM = 20

    libdir = dirname(dirname(realpath(treeshrink_file)))
    tempdir = set_tmp_dir(args["tempdir"])  

    print("Testing R and BMS installation ...")
    if not test_Rlib(libdir):
        print("Failed sanity check on R and BMS installation. Please check your R and BMS version")    
        return

    
    quantiles = [ q for q in args["quantiles"].split()] if args["quantiles"] else ["0.05"]
    
    minImpact = (float(args["minImpact"])/100)+1 if args["minImpact"] else 1.05
    
    scaling = [int(x) for x in args["kscaling"].split(",")] if  args["kscaling"] else [5,2]

    # read gene-to-species mapping file
    g2sp = {}
    if args["g2sp"]:
        print("Reading gene-to-species mapping file")
        with open(args["g2sp"],'r') as fin:
            lineNum = 0
            for line in fin:
                try:
                    g,sp = line.strip().split()
                    g2sp[g] = sp
                except:
                    print("WARNING: failed to parse gene-to-species mapping file on line" + str(lineNum+1))
                    print(line)
                lineNum += 1        

    # exception species
    exceptions = []
    if args["exceptions"]:
        if exists(args["exceptions"]):
            with open(args["exceptions"],'r') as f:
                exceptions = f.read().split()
        else:
            exceptions = args["exceptions"].split()
    exceptions = set(exceptions)

    if args["indir"]:
        treename = splitext(args["tree"])[0]
        subdirs = [d for d in listdir(args["indir"]) if exists(normpath(join(args["indir"],d,args["tree"])))] #if args["tree"] else "input.tre")))]
        #intrees = get_tmp_file(treename + ".trees")
        #with open(intrees,'w') as fout:
        tree_strs = []
        for d in subdirs:
            #treename = args["tree"] if args["tree"] else "input.tre"
            treefile = normpath(join(args["indir"],d,args["tree"]))
            if exists(treefile):
                tree_strs.append(open(treefile,'r').read())
                #fout.write(open(treefile,'r').read())               
        gene_names = [basename(d) for d in subdirs]            
    else:
        #intrees = args["tree"]
        tree_strs = open(args["tree"],'r').readlines()
        gene_names = []

    mode = args["mode"] if args["mode"] else 'auto'
    k = int(args["k"]) if args["k"] else None

    if args["outdir"]:
        outdir = args["outdir"] 
    elif args["indir"]:
        outdir = args["indir"]
    else:
        outdir = splitext(args["tree"])[0] + "_treeshrink"
    if not make_dir(outdir) and args["force"]:
        print("Warning: the output directory " + outdir + " already exists. With --force, all existing files with prefix '" + args["outprefix"]  + "' will be overrided")

    #trees = TreeList.get(path=intrees,schema='newick',preserve_underscores=True)
    #with open(intrees,'r') as f_tree:
    #    tree_strs = f_tree.readlines()
    ntrees = len(tree_strs) 
    if not gene_names:
        gene_names = [str(i) for i in range(ntrees)]

    if mode=='auto' and ntrees < MIN_TREE_NUM:
        print("There are only " + str(ntrees) + " gene trees in the dataset.")
        print("TreeShrink will run in 'All-genes' mode")
        mode='all-genes'

    gene_list = [[] for i in range(ntrees)]
    species_map = {}
    occ = {}
    removing_sets = [ [ [ ] for i in range(ntrees) ] for j in range(len(quantiles)) ]

    for t,a_str in enumerate(tree_strs):
        a_tree = Tree.get(data=a_str,schema='newick',preserve_underscores=True)
        # solve k-shrink
        a_filter = TreeFilter(ddpTree=a_tree,centroid_reroot=args["centroid"],scaling=scaling)
        a_filter.optFilter(d=k)

        # compute species feature (i.e. the max ratio associated with each species for this gene tree)
        mapping = {}
        for i in range(1,len(a_filter.min_diams)):
            if a_filter.min_diams[i] == 0:
                print("Warning: tree %d has no diameter (has only zero branch lengths) after removing %d sequences." %(t+1,i))
                break
            r = a_filter.min_diams[i-1]/a_filter.min_diams[i]
            removals = a_filter.list_removals(d=i)
            for s in removals:
                mapping[s] = r if s not in mapping else max(mapping[s],r)
        # gather per-species distributions and per-gene species features
        for x in mapping:
            s = g2sp[x] if x in g2sp else x
            if mode == 'per-species' or mode == 'auto':
                species_map[s] = [mapping[x]] if s not in species_map else species_map[s]+[mapping[x]]
            #if mode == 'per-species' or mode == 'all-genes' or mode == 'auto':
            gene_list[t].append((x,mapping[x]))
        
        # fit kernel density to this gene's species features (per-gene mode)
        if mode == 'per-gene':
            filename = get_tmp_file("gene_%s.dat" %str(t))
            with open(filename,'w') as f:
                for s in mapping:
                    f.write(str(mapping[s]) + " " +s)
                    f.write("\n")
            if len(mapping) > 1:
                for i,q in enumerate(quantiles):
                    threshold = float(check_output(["Rscript",normpath(join(libdir,"R_scripts","find_threshold_loglnorm.R")),filename,q]).lstrip().rstrip()[4:]) 
                    for s in mapping:
                        if mapping[s] > threshold: 
                            removing_sets[i][t].append(s)
            else:
                print("WARNING: gene " + str(t) + " has too few taxa and is skipped")                
        # update taxon occupancy (only for per-species mode)
        if mode == 'per-species' or mode == 'auto':
            for n in a_tree.leaf_node_iter():
                lb = n.taxon.label
                s = g2sp[lb] if lb in g2sp else lb
                occ[s] = 1 if not s in occ else occ[s]+1
    
    if mode == 'auto' or mode == 'per-species':
        flag = False
        for s in occ:
            if occ[s] < MIN_OCC:
                print ("Species " + s + " only exists in " + str(occ[s]) + " gene trees")
                flag = True
        if flag:
            if mode == 'auto':
                mode = 'all-genes'
                print ("There are species with low occupancy in the dataset. TreeShrink will run in 'All-genes' mode")
            else:
                print ("WARNING: 'Per-species' mode was selected for a dataset having low occupancy species. Consider switching to 'All-genes' mode")
        elif mode == 'auto':
            mode = 'per-species'
            print("TreeShrink will run in 'Per-species' mode ...    ")

# fit kernel density to the per-species distributions and compute per-species threshold (per-species mode)
    if mode == 'per-species':
        for s in sorted(species_map):
            l = len(species_map[s])
            for i in range(occ[s]-l):
                species_map[s].append(1)
            filename = get_tmp_file(s + ".dat")
            with open(filename,'w') as f:
                for v in species_map[s]:
                    f.write(str(v))
                    f.write("\n")
        #if mode == 'per-species':
            thresholds = [ 0 for i in range(len(quantiles)) ]        
            for i,q in enumerate(quantiles):
                thresholds[i] = max(minImpact,float(check_output(["Rscript",normpath(join(libdir,"R_scripts","find_threshold_lkernel.R")),libdir,filename,q]).lstrip().rstrip()[5:]))
                if s not in exceptions:
                    print("%s:\n\t will be cut in %d trees where its impact is above %f for quantile %s" %(s,sum(1 for x in species_map[s] if x>thresholds[i]),thresholds[i],q,))
            species_map[s] = (species_map[s],thresholds)
    #if mode == 'per-species':
        for t,gene in enumerate(gene_list):
            for x,r in gene:
                s = g2sp[x] if x in g2sp else x
                for i,threshold in enumerate(species_map[s][1]):
                    if r > threshold:
                        removing_sets[i][t].append(x)
                    

# fit kernel density to all the species features across all genes and compute the global threshold (all-gene mode) 
    if mode == 'all-genes':
        filename = get_tmp_file("all_genes" + ".dat")
        with open(filename,'w') as f:
            for gene in gene_list:
                for s,r in gene:
                    f.write(str(r))
                    f.write("\n")
        for i,q in enumerate(quantiles):
            threshold = float(check_output(["Rscript",normpath(join(libdir,"R_scripts","find_threshold_lkernel.R")),libdir,filename,q]).lstrip().rstrip()[5:])
            for t,gene in enumerate(gene_list):
                for x,r in gene:
                    #s = g2sp[x] if x in g2sp else x
                    if r > threshold:
                        removing_sets[i][t].append(x)

    print("Writing output ...\n")

    fName,ext = splitext(basename(args["tree"]))
    ext = ext if ext else '.nwk'
    prefix = args["outprefix"]
    counter = 0
    # check if the outdir or any of its subdirs already has files with the specified prefix
    for File in listdir(outdir):
        if File.startswith(prefix):
             search_counter = re.search(r'\d+', File[len(prefix):])
             counter = max(counter,1 if not search_counter else int(search_counter.group())+1)
        if isdir(normpath(join(outdir,File))):
            for File1 in listdir(normpath(join(outdir,File))):
                if File1.startswith(prefix):
                    search_counter = re.search(r'\d+', File1[len(prefix):])
                    counter = max(counter,1 if not search_counter else int(search_counter.group())+1)
    if counter >0 and not args["force"]:
        print("WARNING: " + outdir + " has already had some files with prefix '" + prefix + "'. Automatically changes prefix to '" + prefix + str(counter) + "' to avoid overriding. Rerun with --force if you wish to override existing files.")            
        prefix = prefix + str(counter)

    # write summary file
    filename= normpath(join(outdir,prefix + "_summary.txt"))                
    with open(filename,'w') as f:
        f.write("Gene Species Taxon Signature\n")
        for t,gene in enumerate(gene_list):
            for x,r in gene:
                s = g2sp[x] if x in g2sp else x
                f.write(gene_names[t] + " " + s + " " + x + " " + str(log(r)))
                f.write("\n")
    
    # write shrunk trees, removing sets, and filtered alignments
    # Dendropy's filter_leaf_nodes() seems to have problem
    # i.e. it produces the trees that the treecmp tool cannot compute the MS distance (need further exploration)
    # use home-made code to prune the tree instead
     
    for i,RS in enumerate(removing_sets):
        #trees_shrunk = deepcopy(trees)
        RS_tag = '' if (len(removing_sets) < 2) else '_' + quantiles[i]
        tree_tag = '' if (len(removing_sets) < 2) else '_' + quantiles[i]
        aln_tag = '' if (len(removing_sets) < 2) else '_' + quantiles[i]
        
        if args["indir"] is None:
            outfile = normpath(join(outdir,prefix + RS_tag + ".txt"))
            with open(outfile,'w') as f:
                for item in RS:
                    for s in item:
                        if s not in exceptions:
                            f.write(s + "\t")
                    f.write("\n")
            for a_str,rs in zip(tree_strs,RS):
                tree = Tree.get(data=a_str,schema='newick',preserve_underscores=True)
                prune_tree(tree,set(rs)-exceptions)
                tree_as_newick(tree,outfile=normpath(join(outdir,prefix + tree_tag + ext)),append=True)
        else:
            for sd,item in zip(subdirs,RS):
                make_dir(normpath(join(outdir,sd)))
                outfile = normpath(join(outdir,sd, prefix + RS_tag + ".txt"))
                with open(outfile,'w') as f:
                    for s in item:
                        if s not in exceptions:
                            f.write(s + "\t")
            for sd,a_str,rs in zip(subdirs,tree_strs,RS):
                tree = Tree.get(data=a_str,schema='newick',preserve_underscores=True)
                L = set(x.taxon.label for x in tree.leaf_node_iter())
                rs1 = set(rs)-exceptions
                prune_tree(tree,rs1)
                treefile = normpath(join(outdir,sd, prefix + tree_tag + ext))
                #tree.write_to_path(treefile,'newick',unquoted_underscores=True,real_value_format_specifier=".16g")
                tree_as_newick(tree,outfile=treefile,append=False)
                
                aln_filename = args["alignment"] if args["alignment"] else "input.fasta"
                alnName,alnExt = splitext(aln_filename)
                alnExt = alnExt if alnExt else '.fasta'
                input_aln = normpath(join(args["indir"],sd,aln_filename))
                if isfile(input_aln): 
                    output_aln = normpath(join(outdir,sd,prefix+aln_tag+alnExt))
                    alg = CompactAlignment()
                    alg.read_file_object(input_aln,'fasta')
                    S=set(alg.keys())
                    if (L.difference(alg.keys())) or S.difference(L):
                        print("ERROR: For gene %s, alignment names don't match tree names. Will skip it.\n\tonly in tree:\t%s\n\tonly in alignment:\t%s"%(sd,str(L.difference(S)),str(S.difference(L))))
                    else:
                        alg.remove_all(rs1)
                        alg.mask_gapy_sites(1)
                        alg.write(output_aln,'fasta')
    if not args["tempdir"]:
        rmtree(tempdir)

    print("Output files written to " + outdir + " with prefix " + prefix + ".") 

    
if __name__ == "__main__":
    main()
