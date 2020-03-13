#!/usr/bin/env python2

import sys
import numpy as np
import scipy.spatial.distance as spd 
import scipy.cluster.hierarchy as sph
from scipy import stats
import matplotlib
#matplotlib.use('Agg')
import pylab
import pandas as pd
import matplotlib.pyplot as plt

class ReadCmd:

    def __init__( self ):
        import argparse as ap
        import textwrap

        p = ap.ArgumentParser( description= "TBA" )
        arg = p.add_argument
        
        arg( '-i', '--inp', '--in', metavar='INPUT_FILE', type=str, nargs='?', default=sys.stdin,
             help= "The input matrix" )
        arg( '-o', '--out', metavar='OUTPUT_FILE', type=str, nargs='?', default=None,
             help= "The output image file [image on screen of not specified]" )

        arg( '-m', '--metadata_file', type=str, default='None',
             help= "The input metadata file [default None]" )

        DataMatrix.input_parameters( p )
        BarPlot.input_parameters( p )
        self.args  = p.parse_args()

    def check_consistency( self ):
        pass

    def get_args( self ):
        return self.args

class DataMatrix:
    datatype = 'data_matrix'
    
    @staticmethod
    def input_parameters( parser ):
        dm_param = parser.add_argument_group('Input data matrix parameters')
        arg = dm_param.add_argument

        arg( '--sep', type=str, default='\t' )
        arg( '-f', '--feat', type=str, default=None, required = True,
             help = "Name of the feature to plot"
                    "[or the ending string if --endswith is specified]")
        arg( '--endswith', action='store_true',
             help = "Match the ending part of the feature name" )
        arg( '--fname_row', type=int, default=0,
             help = "row number containing the names of the features "
                    "[default 0, specify -1 if no names are present in the matrix")
        arg( '--sname_row', type=int, default=0,
             help = "column number containing the names of the samples "
                    "[default 0, specify -1 if no names are present in the matrix")
        arg( '--skip_rows', type=str, default=None,
             help = "Row numbers to skip (0-indexed, comma separated) from the input file"
                    "[default None, meaning no rows skipped")
        arg( '--def_na', type=float, default=None,
             help = "Set the default value for missing values [default None which means no replacement]")

    def __init__( self, input_file, args ):
        self.args = args
        toskip = [int(l) for l in self.args.skip_rows.split(",")]  if self.args.skip_rows else None
        self.table = pd.read_table( 
                input_file, sep = self.args.sep, skipinitialspace = True, skiprows = toskip,
                                  header = self.args.fname_row if self.args.fname_row > -1 else None,
                                  index_col = self.args.sname_row if self.args.sname_row > -1 else None
                                    )

        rows = []

        if self.args.endswith:
            for n in self.table.index:
                if n.endswith( self.args.feat  ):
                    rows.append( n )
        elif self.args.feat in self.table.index:
            rows.append( self.args.feat )
        self.table = self.table.reindex( index=rows )

        if not len(rows):
            sys.stderr.write("Error, feat "+self.args.feat+" not found!")
            sys.exit()
        if len(rows) > 1:
            sys.stderr.write("Error, multiple features matching "+self.args.feat+" !")
            sys.exit()

        if not self.args.def_na is None:
            self.table = self.table.fillna( self.args.def_na )

    def get_numpy_matrix( self ): 
        return self.table
    
    def get_snames( self ):
        return list(self.table.index)
    
    def get_fnames( self ):
        return list(self.table.columns)
   
    def save_matrix( self, output_file ):
        self.table.to_csv( output_file, sep = '\t' )

class MetadataMatrix:
    datatype = 'metadata_matrix'
    
    @staticmethod
    def input_parameters( parser ):
        dm_param = parser.add_argument_group('Input metadata file')
        arg = dm_param.add_argument

        arg( '--sep', type=str, default='\t' )
        arg( '--fname_row', type=int, default=0,
             help = "row number containing the names of the features "
                    "[default 0, specify -1 if no names are present in the matrix")
        arg( '--def_na', type=float, default=None,
             help = "Set the default value for missing values [default None which means no replacement]")

    def __init__( self, input_file, args ):
        self.args = args
        self.table = pd.read_table( 
                input_file, sep = self.args.sep, skipinitialspace = True, 
                        #header = self.args.fname_row if self.args.fname_row > -1 else None,
                                  index_col = self.args.sname_row if self.args.sname_row > -1 else None
                                    )

        if not self.args.def_na is None:
            self.table = self.table.fillna( self.args.def_na )
    
    def get_snames( self ):
        return list(self.table.index)
    
    def get_fnames( self ):
        return list(self.table.columns)

    def get_table( self ):
        return self.table

class BarPlot:
    datatype = 'barplot'

    @staticmethod
    def input_parameters( parser ):
        hm_param = parser.add_argument_group('Heatmap options')
        arg = hm_param.add_argument

        arg( '--dpi', type=int, default=72,
             help = "Image resolution in dpi [default 72]")
        arg( '-C', '--color_condition', type=str, default=None,
             help = "The name of the metadata column used for coloring")
        arg( '-H', '--hatch_condition', type=str, default=None,
             help = "The name of the metadata column used for hatching")
        arg( '-G', '--group_condition', type=str, default=None,
             help = "The name of the metadata column used for grouping")
        arg( '-t', '--title', type=str, default=None,
             help = "The title of the plot [default no title]")
        arg( '-l', '--log_scale', action='store_true',
             help = "Log scale" )

    
    def __init__( self, numpy_matrix, metadata_matrix, args = None ):
        self.numpy_matrix = numpy_matrix
        self.mmatrix = metadata_matrix
        self.args = args

    def draw( self ):

        fig = plt.figure( figsize=(20,8)  )
        ax = fig.add_subplot(111)

        width = 0.65      

        names = list(self.numpy_matrix.index)
        n0 = names[0]

        tp = self.numpy_matrix.to_dict()
        
        keys = sorted(tp)
        
        if self.args.color_condition not in self.mmatrix:
            self.args.color_condition = None
        cond_values = [None] if self.args.color_condition is None else sorted(set(self.mmatrix[self.args.color_condition]) )
        if self.args.hatch_condition not in self.mmatrix:
            self.args.hatch_condition = None
        hatch_values = [None] if self.args.hatch_condition is None else sorted(set(self.mmatrix[self.args.hatch_condition]) )
       
        if self.args.group_condition:
            group_values = list(sorted(set(self.mmatrix[self.args.group_condition])))
            keys = sorted( keys, key=lambda x:group_values.index(self.mmatrix[self.args.group_condition][x]) )
        else:
            keys, group_values = sorted( keys ), []

        ind = np.arange( len(tp) )
        pos = ind-width/2

        hatches = ['//','\\\\','++','--','xx']
        cols = ['r','g','c','b']
        minv,maxv = 0.0, max([v[n0] for v in tp.values()])
        
        bar_sets = []
        for i,c in (enumerate(cond_values) if len(cond_values) > 0 else None): 
            for j,h in enumerate(hatch_values): 
                values = [(tp[k][n0] if (c is None or self.mmatrix[self.args.color_condition][k] == c) 
                                    and (h is None or self.mmatrix[self.args.hatch_condition][k] == h) else 0.0) for k in keys]
                b = ax.bar(pos, values, width, hatch=hatches[j%len(hatches)] if len(hatch_values) > 1 else "", color=cols[i%len(cols)])
                cond = self.args.color_condition + " "+str(c).strip()+", " if c else ""
                hatch = self.args.hatch_condition + " "+str(h).strip()+", " if h else ""
                bar_sets.append( (b,cond+hatch) )

        v0 = ind[0]-0.5
        vm1 = v0
        ax.plot([v0,v0],[minv,maxv],"--",linewidth=2,color='k')
        for g in group_values:
            vm1 = v0
            v0 += list(self.mmatrix[self.args.group_condition]).count(g)
            ax.plot([v0,v0],[minv,maxv],"--",linewidth=2,color='k')
            ax.text( (vm1+v0)*0.5, maxv * 0.9, str(g), horizontalalignment='center', verticalalignment='center' )
            #ax.text( (vm1+v0)*0.5, maxv * 0.9, str(round(g,1)), horizontalalignment='center', verticalalignment='center' )

        if self.args.color_condition or self.args.hatch_condition:
            leg = ax.legend( zip(*bar_sets)[0], zip(*bar_sets)[1], bbox_to_anchor=(1.02, 0,0.3,1), loc=1,
                           ncol=1, mode="expand", borderaxespad=0., frameon = False)

        ax.set_xlim(-width,ind[-1]+width)
        ax.set_ylim(0,maxv)
        ax.set_xticks( ind )
        ax.set_xticklabels( keys, rotation = 90 )
        ax.set_title( self.args.title or "" )

        if not self.args.out:
            plt.show()
        else:
            fig.savefig( self.args.out, bbox_inches='tight', dpi = self.args.dpi, 
                         bbox_extra_artists=((fig.get_axes()[0].get_legend(),) if self.args.color_condition or self.args.hatch_condition else None) ) #dpi = self.args.dpi )

if __name__ == '__main__':

    read = ReadCmd( )
    read.check_consistency()
    args = read.get_args()

    dm = DataMatrix( args.inp, args )
    mdm = MetadataMatrix( args.metadata_file, args ) 

    bp = BarPlot( dm.get_numpy_matrix(), mdm.get_table(),args )      
    bp.draw()




