#!/usr/bin/env python

import argparse as ap
import sys
import math

__author__ = 'Nicola Segata (nicola.segata@unitn.it)'
__version__ = '0.9'
__date__ = '17th March 2013'

def read_params(args):
    p = ap.ArgumentParser( description= 
                           
            "DESCRIPTION\n"
            "  metahplan2graphlan.py version "+__version__+" ("+__date__+")\n" 
            "  Convert MetaPhlAn outputs to GraPhlAn imput format. \n"
            "  AUTHORS: "+__author__+"\n\n"
            "EXAMPLE\n"
            "  metaphlan2graphlan.py metahplan-out.txt --tree_file out_graphlan_tree.txt --annot_file out_annotation_options.txt",
            formatter_class=ap.RawTextHelpFormatter )
    arg = p.add_argument
    arg( 'inp', metavar='INPUT_FILE', type=str, nargs='?', default=None, help= 
         "Merged MetaPhlAn data.\n")   
    arg( '--tree_file', type=str, required = True, help = "The output tree file for GraPhlAn" ) 
    arg( '--annot_file', type=str, required = True, help = "The annotation file for GraPhlAn" ) 
    arg( '--max_annot_clades', type=int, default = 10, help = "The maximum number of clades to annotate [default 10]" ) 
    arg( '--min_annot_lev', type=int, default = 2, help = "The minimum number of levels required to label a clade [default 2, meaning microbial orders]" ) 
    arg( '--max_annot_lev', type=int, default = 5, help = "The maximum number of levels required to label a clade [default 5, meaning microbial species]" ) 
    arg( '--ext_keys_start_lev', type=int, default = 5, help = "The level at which annotations are added using external legend keys [default 5, meaning microbial species]" ) 
    arg( '--coloring_lev', type=int, default = 3, help = "The level used for color differentiation [default 3, meaning microbial families]" ) 
    return vars(p.parse_args())


def readMetaPhlAnFile(fileName):
    observation_ids = []
    observation_md = []
    data = {} 
    sample_md = None
    
    f = open(fileName, "r")
    fileData = f.read()
    f.close()
    
    lines = fileData.split("\n")

    header = lines[0].split("\t")
        
    headerLen = len(header)
    sample_ids = header[1:]
    datalines = lines
    
    for i in range(len(datalines)):
        line = datalines[i].split("\t")
        if line[0] == "ID":
            continue
        if (datalines[i] == ""): continue
        if len(line) != headerLen: 
            raise Exception("readMetaPhlAnFile(): Unexpected input file format (lines do not match header). Is the input file generated with merge_metaphlan_tables.py?")

        taxonomy = line[0].split("|")
        datarow = line[1:]
        
        datarowScaled = []
        for v in datarow:
            datarowScaled.append(float(v))

        data[line[0].replace("|",".")] = 10.0+100.0*math.log(1.0+sum(datarowScaled)/float(len(datarowScaled)))

    return data

cols = ["#0000cd","#006400","#b22222","#00bfff","#00ff7f","#ffd700","#8a2be2"]

if __name__ == '__main__':
    pars = read_params( sys.argv )
    
    try:
        if pars['inp'] is None:
            raise Exception("Please give MetaPhlAn file (merged) to be converted.")

        taxa2avg_ab = readMetaPhlAnFile(pars['inp'])

    except Exception as exception:
        sys.stderr.write("Error: {0}\n".format(str(exception)))

    with open(pars['tree_file'],"w") as outf:
        outf.write( "\n".join( taxa2avg_ab ) )

    col_levs = []
    with open(pars['annot_file'],"w") as outf:
        i = 0
        for t,a in sorted(taxa2avg_ab.items(),key=lambda x:x[1],reverse=True):
            levs = t.count(".")
            outf.write( "\t".join( [t,"clade_marker_size",str(a)] ) +"\n" )
            if i < pars['max_annot_clades'] and levs > pars['min_annot_lev'] and levs <= pars['max_annot_lev']:
                i += 1
                annot = "*" if levs < pars['ext_keys_start_lev'] else "*:*"
                outf.write( "\t".join( [t,"annotation",annot] ) +"\n" )
                cl = t.split(".")[pars['coloring_lev']]
                if cl not in col_levs:
                    col_levs.append(cl)
                col = cols[col_levs.index( cl ) % len(cols)]
                outf.write( "\t".join( [t,"clade_marker_color",col] ) +"\n" ) 
                outf.write( "\t".join( [t,"annotation_background_color",col] ) +"\n" ) 







