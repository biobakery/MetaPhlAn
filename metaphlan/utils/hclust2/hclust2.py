#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.ticker as ticker
import scipy.spatial.distance as spd
import scipy.cluster.hierarchy as sph
from scipy import stats
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.fonttype'] = 'none'
import pylab
import pandas as pd
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
try:
    import cPickle as pickle
except:
    import pickle # Python3 compatible (Lauren McIver)
import math

sys.setrecursionlimit(10000)

# samples on rows

class SqrtNorm(matplotlib.colors.Normalize):
    """
    Normalize a given value to the 0-1 range on a square root scale
    """
    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        result = np.ma.masked_less_equal(result, 0, copy=False)

        self.autoscale_None(result)
        vmin, vmax = self.vmin, self.vmax
        if vmin > vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        elif vmin <= 0:
            raise ValueError("values must all be positive")
        elif vmin == vmax:
            result.fill(0)
        else:
            if clip:
                mask = np.ma.getmask(result)
                result = np.ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)
            # in-place equivalent of above can be much faster
            resdat = result.data
            mask = result.mask
            if mask is np.ma.nomask:
                mask = (resdat <= 0)
            else:
                mask |= resdat <= 0
            np.copyto(resdat, 1, where=mask)
            np.sqrt(resdat, resdat)
            resdat -= np.sqrt(vmin)
            resdat /= (np.sqrt(vmax) - np.sqrt(vmin))
            result = np.ma.array(resdat, mask=mask, copy=False)
        if is_scalar:
            result = result[0]
        return result

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        vmin, vmax = self.vmin, self.vmax

        if matplotlib.cbook.iterable(value):
            val = np.ma.asarray(value)
            return vmin * np.ma.power((vmax / vmin), val)
        else:
            return vmin * pow((vmax / vmin), value)

    def autoscale(self, A):
        '''
        Set *vmin*, *vmax* to min, max of *A*.
        '''
        A = np.ma.masked_less_equal(A, 0, copy=False)
        self.vmin = np.ma.min(A)
        self.vmax = np.ma.max(A)

    def autoscale_None(self, A):
        ' autoscale only None-valued vmin or vmax'
        if self.vmin is not None and self.vmax is not None:
            return
        A = np.ma.masked_less_equal(A, 0, copy=False)
        if self.vmin is None:
            self.vmin = np.ma.min(A)
        if self.vmax is None:
            self.vmax = np.ma.max(A)

class DataMatrix:
    datatype = 'data_matrix'

    @staticmethod
    def input_parameters( parser ):
        dm_param = parser.add_argument_group('Input data matrix parameters')
        arg = dm_param.add_argument

        arg( '--sep', type=str, default='\t' )
        arg( '--out_table', type=str, default=None,
             help = 'Write processed data matrix to file' )
        arg( '--fname_row', type=int, default=0,
             help = "row number containing the names of the features "
                    "[default 0, specify -1 if no names are present in the matrix")
        arg( '--sname_row', type=int, default=0,
             help = "column number containing the names of the samples "
                    "[default 0, specify -1 if no names are present in the matrix")
        arg( '--metadata_rows', type=str, default=None,
             help = "Row numbers to use as metadata"
                    "[default None, meaning no metadata")
        arg( '--skip_rows', type=str, default=None,
             help = "Row numbers to skip (0-indexed, comma separated) from the input file"
                    "[default None, meaning no rows skipped")
        arg( '--sperc', type=int, default=90,
             help = "Percentile of sample value distribution for sample selection" )
        arg( '--fperc', type=int, default=90,
             help = "Percentile of feature value distribution for sample selection" )
        arg( '--stop', type=int, default=None,
             help = "Number of top samples to select (ordering based on percentile specified by --sperc)" )
        arg( '--ftop', type=int, default=None,
             help = "Number of top features to select (ordering based on percentile specified by --fperc)" )
        arg( '--def_na', type=float, default=None,
             help = "Set the default value for missing values [default None which means no replacement]")

    def __init__( self, input_file, args ):
        self.args = args
        self.metadata_rows =  []
        self.metadata_table = None
        toskip = [int(l) for l in self.args.skip_rows.split(",")]  if self.args.skip_rows else []
        if self.args.metadata_rows:
            self.metadata_rows = list([int(a) for a in self.args.metadata_rows.split(",")])
            mdr = self.metadata_rows[::]
            for t in toskip:
                for i,m in enumerate(mdr):
                    if t <= m:
                        self.metadata_rows[i] -= 1
        if self.metadata_rows:
            header = [self.args.fname_row]+self.metadata_rows if self.args.fname_row > -1 else self.metadata_rows
        else:
            header = self.args.fname_row if self.args.fname_row > -1 else None
        self.table = pd.read_table(
                input_file, sep = self.args.sep, # skipinitialspace = True,
                                  skiprows = sorted(toskip) if isinstance(toskip, list) else toskip,
                                  header = sorted(header) if isinstance(header, list) else header,
                                  index_col = self.args.sname_row if self.args.sname_row > -1 else None
                                    )

        def select( perc, top  ):
            self.table['perc'] = self.table.apply(lambda x: stats.scoreatpercentile(x,perc),axis=1)

            if top <= len(self.table['perc']):
                m = sorted(self.table['perc'])[-top]
            else:
                print('W ftop param value (' + str(top) + ') out of bound (len:' + str(len(self.table['perc'])) + '). Selecting all the values from input.')
                m = sorted(self.table['perc'])[0]

            self.table = self.table[self.table['perc'] >= m ]
            del self.table['perc']

        if not self.args.def_na is None:
            self.table = self.table.fillna( self.args.def_na )

        if self.args.ftop:
            select( self.args.fperc, self.args.ftop )

        if self.args.stop:
            self.table = self.table.T
            select( self.args.sperc, self.args.stop )
            self.table = self.table.T


        # add missing values

    def get_numpy_matrix( self ):
        return np.matrix(self.table)

    #def get_metadata_matrix( self ):
    #    return self.table.columns

    def get_snames( self ):
        #return list(self.table.index)
        return self.table.columns

    def get_fnames( self ):
        #print self.table.columns.names
        #print self.table.columns
        return list(self.table.index)

    def get_averages(self, by_row = True) :
        return self.table.mean(axis = 1 if by_row else 0)

    def save_matrix( self, output_file ):
        self.table.to_csv( output_file, sep = '\t' )

class DistMatrix:
    datatype = 'distance_matrix'

    @staticmethod
    def input_parameters( parser ):
        dm_param = parser.add_argument_group('Distance parameters')
        arg = dm_param.add_argument

        dist_funcs = [  "euclidean","minkowski","cityblock","seuclidean",
                        "sqeuclidean","cosine","correlation","hamming",
                        "jaccard","chebyshev","canberra","braycurtis",
                        "mahalanobis","yule","matching","dice",
                        "kulsinski","rogerstanimoto","russellrao","sokalmichener",
                        "sokalsneath","wminkowski","ward" ]

        arg( '--f_dist_f', type=str, default="correlation",
             help = "Distance function for features [default correlation]")
        arg( '--s_dist_f', type=str, default="euclidean",
             help = "Distance function for sample [default euclidean]")
        arg( '--load_dist_matrix_f', type=str, default=None,
             help = "Load the distance matrix to be used for features [default None].")
        arg( '--load_dist_matrix_s', type=str, default=None,
             help = "Load the distance matrix to be used for samples [default None].")
        arg( '--load_pickled_dist_matrix_f', type=str, default=None,
             help = "Load the distance matrix to be used for features as previously saved as pickle file using hclust2 itself [default None].")
        arg( '--load_pickled_dist_matrix_s', type=str, default=None,
             help = "Load the distance matrix to be used for samples as previously saved as pickle file using hclust2 itself [default None].")
        arg( '--save_pickled_dist_matrix_f', type=str, default=None,
             help = "Save the distance matrix for features to file [default None].")
        arg( '--save_pickled_dist_matrix_s', type=str, default=None,
             help = "Save the distance matrix for samples to file [default None].")

    def __init__( self, data, args = None ):
        self.sdf = args.s_dist_f
        self.fdf = args.f_dist_f

        self.s_cdist_matrix, self.f_cdist_matrix = None, None

        self.numpy_full_matrix = (data if
                type(data) == np.matrixlib.defmatrix.matrix else None)

    def compute_f_dists( self ):
        if args.load_pickled_dist_matrix_f:
            with open( args.load_pickled_dist_matrix_f ) as inp:
                self.f_cdist_matrix = pickle.load( inp )
        elif args.load_dist_matrix_f:
            self.f_cdist_matrix = spd.squareform( np.matrix( pd.read_table( args.load_dist_matrix_f, sep ='\t', index_col = None, header = None  ) ) )
        else:
            dt = self.numpy_full_matrix

            if self.fdf == "spearman":
                dt_ranked = np.matrix([stats.rankdata(d) for d in dt])
                self.f_cdist_matrix = spd.pdist( dt_ranked, "correlation" )
                return

            if self.fdf == 'mhamming':
                dt_ranked = np.matrix([[(0 if l == 0 else 1) for l in np.nditer(d)] for d in dt])
                self.f_cdist_matrix = spd.pdist( dt_ranked, "hamming" )
                return

            if self.fdf == 'lbraycurtis':
                dt_ranked = np.matrix([[(math.log(l) if l else 0.0) for l in np.nditer(d)] for d in dt])
                self.f_cdist_matrix = spd.pdist( dt_ranked, "braycurtis" )
                return

            if self.fdf == "pearson":
                self.fdf = 'correlation'

            self.f_cdist_matrix = spd.pdist( dt, self.fdf )

        if args.save_pickled_dist_matrix_f:
            with open( args.save_pickled_dist_matrix_f, "wb" ) as outf:
                pickle.dump( self.f_cdist_matrix, outf )

    def compute_s_dists( self ):
        if args.load_pickled_dist_matrix_s:
            with open( args.load_pickled_dist_matrix_s ) as inp:
                self.s_cdist_matrix = pickle.load( inp )
        elif args.load_dist_matrix_s:
            self.s_cdist_matrix = spd.squareform( np.matrix( pd.read_table( args.load_dist_matrix_s, sep ='\t', index_col = None, header = None  ) ) )
        else:
            done = False
            dt = self.numpy_full_matrix.transpose()

            if self.sdf == "spearman":
                dt_ranked = np.matrix([stats.rankdata(d) for d in dt])
                self.s_cdist_matrix = spd.pdist( dt_ranked, "correlation" )
                done = True

            if self.sdf == 'mhamming':
                dt_ranked = np.matrix([[(0 if l == 0 else 1) for l in np.nditer(d)] for d in dt])
                self.s_cdist_matrix = spd.pdist( dt_ranked, "hamming" )
                done = True

            if self.sdf == 'lbraycurtis':
                dt_ranked = np.matrix([[(math.log(l) if l else 0.0) for l in np.nditer(d)] for d in dt])
                self.s_cdist_matrix = spd.pdist( dt_ranked, "braycurtis" )
                done = True

            if self.sdf == 'sbraycurtis':
                dt_ranked = np.matrix([[(math.sqrt(l) if l else 0.0) for l in np.nditer(d)] for d in dt])
                self.s_cdist_matrix = spd.pdist( dt_ranked, "braycurtis" )
                done = True

            if self.sdf == "pearson":
                self.sdf = 'correlation'

            if not done:
                self.s_cdist_matrix = spd.pdist( dt, self.sdf )

        if args.save_pickled_dist_matrix_s:
            with open( args.save_pickled_dist_matrix_s, "wb" ) as outf:
                pickle.dump( self.s_cdist_matrix, outf )

    def get_s_dm( self ):
        return self.s_cdist_matrix

    def get_f_dm( self ):
        return self.f_cdist_matrix

class HClustering:
    datatype = 'hclustering'

    @staticmethod
    def input_parameters( parser ):
        cl_param = parser.add_argument_group('Clustering parameters')
        arg = cl_param.add_argument

        linkage_method = [ "single","complete","average",
                           "weighted","centroid","median",
                           "ward" ]
        arg( '--no_fclustering', action='store_true',
             help = "avoid clustering features" )
        arg( '--no_plot_fclustering', action='store_true',
             help = "avoid plotting the feature dendrogram" )
        arg( '--no_sclustering', action='store_true',
             help = "avoid clustering samples" )
        arg( '--no_plot_sclustering', action='store_true',
             help = "avoid plotting the sample dendrogram" )
        arg( '--flinkage', type=str, default="average",
             help = "Linkage method for feature clustering [default average]")
        arg( '--slinkage', type=str, default="average",
             help = "Linkage method for sample clustering [default average]")

    def get_reordered_matrix( self, matrix, sclustering = True, fclustering = True ):
        if not sclustering and not fclustering:
            return matrix

        idx1 = self.sdendrogram['leaves'] if sclustering else None   # !!!!!!!!!!!
        idx2 = self.fdendrogram['leaves'][::-1] if fclustering else None

        if sclustering and fclustering:
            return matrix[idx2,:][:,idx1]
        if fclustering:
            return matrix[idx2,:][:]
        if sclustering: # !!!!!!!!!!!!
            return matrix[:][:,idx1]

    def get_reordered_sample_labels( self, slabels ):
        return [slabels[i] for i in self.sdendrogram['leaves']]

    def get_reordered_feature_labels( self, flabels ):
        return [flabels[i] for i in self.fdendrogram['leaves']]

    def __init__( self, s_dm, f_dm, args = None ):
        self.s_dm = s_dm
        self.f_dm = f_dm
        self.args = args
        self.sclusters = None
        self.fclusters = None
        self.sdendrogram = None
        self.fdendrogram = None

    def shcluster( self, dendrogram = True ):
        self.shclusters = sph.linkage(self.s_dm, method=args.slinkage)
        if dendrogram:
            self.sdendrogram = sph.dendrogram( self.shclusters, no_plot=True )

    def fhcluster( self, dendrogram = True ):
        self.f_dm = [abs(round(i,15)) for i in self.f_dm]
        self.fhclusters = sph.linkage(self.f_dm, method=args.flinkage)
        if dendrogram:
            self.fdendrogram = sph.dendrogram( self.fhclusters, no_plot=True )

    def get_shclusters( self ):
        return self.shclusters

    def get_fhclusters( self ):
        return self.fhclusters

    def get_sdendrogram( self ):
        return self.sdendrogram

    def get_fdendrogram( self ):
        return self.fdendrogram


class Heatmap:
    datatype = 'heatmap'

    bbcyr = {'red':  (  (0.0, 0.0, 0.0),
                        (0.25, 0.0, 0.0),
                        (0.50, 0.0, 0.0),
                        (0.75, 1.0, 1.0),
                        (1.0, 1.0, 1.0)),
             'green': ( (0.0, 0.0, 0.0),
                        (0.25, 0.0, 0.0),
                        (0.50, 1.0, 1.0),
                        (0.75, 1.0, 1.0),
                        (1.0, 0.0, 1.0)),
             'blue': (  (0.0, 0.0, 0.0),
                        (0.25, 1.0, 1.0),
                        (0.50, 1.0, 1.0),
                        (0.75, 0.0, 0.0),
                        (1.0, 0.0, 1.0))}

    bbcry = {'red':  (  (0.0, 0.0, 0.0),
                        (0.25, 0.0, 0.0),
                        (0.50, 0.0, 0.0),
                        (0.75, 1.0, 1.0),
                        (1.0, 1.0, 1.0)),
             'green': ( (0.0, 0.0, 0.0),
                        (0.25, 0.0, 0.0),
                        (0.50, 1.0, 1.0),
                        (0.75, 0.0, 0.0),
                        (1.0, 1.0, 1.0)),
             'blue': (  (0.0, 0.0, 0.0),
                        (0.25, 1.0, 1.0),
                        (0.50, 1.0, 1.0),
                        (0.75, 0.0, 0.0),
                        (1.0, 0.0, 1.0))}

    bcry = {'red':  (   (0.0, 0.0, 0.0),
                        (0.33, 0.0, 0.0),
                        (0.66, 1.0, 1.0),
                        (1.0, 1.0, 1.0)),
             'green': ( (0.0, 0.0, 0.0),
                        (0.33, 1.0, 1.0),
                        (0.66, 0.0, 0.0),
                        (1.0, 1.0, 1.0)),
             'blue': (  (0.0, 1.0, 1.0),
                        (0.33, 1.0, 1.0),
                        (0.66, 0.0, 0.0),
                        (1.0, 0.0, 1.0))}


    my_colormaps = [    ('bbcyr',bbcyr),
                        ('bbcry',bbcry),
                        ('bcry',bcry)]

    #dcols = ['#ca0000','#0087ff','#00ba1d','#cf00ff','#00dbe2','#ffaf00','#0017f4','#006012','#e175ff','#877878','#050505','#b5cf00','#ff8a8a','#aa6400','#50008a','#00ff58']
    dcols = ['#ca0000','#0087ff','#00ba1d','#cf00ff','#00dbe2','#ffaf00','#0017f4','#006012','#e175ff','#877878','#505050','#b5cf00','#ff8a8a','#aa6400','#50008a','#00ff58','#6F1A1A','#FFCC99','#33FF33','#009999','#CC0066','#99004c','#C0C0C0',"#666600","#CCFF99","#660066","#9370DB","#D8BFD8","#BC8F8F","#2F4F4F","#FF6347","#CD5C5C","#FF0000","#00FF00","#000080"]


    @staticmethod
    def input_parameters( parser ):
        hm_param = parser.add_argument_group('Heatmap options')
        arg = hm_param.add_argument

        arg( '--dpi', type=int, default=150,
             help = "Image resolution in dpi [default 150]")
        arg( '-l', '--log_scale', action='store_true',
             help = "Log scale" )
        arg( '--title', type=str, default=None,
             help = "Title of the plot" )
        arg( '--title_fontsize', type=int, default=10,
             help = "Font size of the title" )
        arg( '-s', '--sqrt_scale', action='store_true',
             help = "Square root scale" )
        arg( '--no_slabels', action='store_true',
             help = "Do not show sample labels" )
        arg( '--minv', type=float, default=None,
             help = "Minimum value to display in the color map [default None meaning automatic]" )
        arg( '--maxv', type=float, default=None,
             help = "Maximum value to display in the color map [default None meaning automatic]" )
        arg( '--no_flabels', action='store_true',
             help = "Do not show feature labels" )
        arg( '--max_slabel_len', type=int, default=25,
             help = "Max number of chars to report for sample labels [default 15]" )
        arg( '--max_flabel_len', type=int, default=25,
             help = "Max number of chars to report for feature labels [default 15]" )
        arg( '--flabel_size', type=int, default=10,
             help = "Feature label font size [default 10]" )
        arg( '--slabel_size', type=int, default=10,
             help = "Sample label font size [default 10]" )
        arg( '--fdend_width', type=float, default=1.0,
             help = "Width of the feature dendrogram [default 1 meaning 100%% of default heatmap width]")
        arg( '--sdend_height', type=float, default=1.0,
             help = "Height of the sample dendrogram [default 1 meaning 100%% of default heatmap height]")
        arg( '--metadata_height', type=float, default=.05,
             help = "Height of the metadata panel [default 0.05 meaning 5%% of default heatmap height]")
        arg( '--metadata_separation', type=float, default=.01,
             help = "Distance between the metadata and data panels. [default 0.001 meaning 0.1%% of default heatmap height]")
        arg( '--colorbar_font_size', type=float, default=12,
             help = "Color bar label font size [default 12]")
        arg( '--image_size', type=float, default=8,
             help = "Size of the largest between width and eight size for the image in inches [default 8]")
        arg( '--cell_aspect_ratio', type=float, default=1.0,
             help = "Aspect ratio between width and height for the cells of the heatmap [default 1.0]")
        col_maps = ['Accent', 'Blues', 'BrBG', 'BuGn', 'BuPu', 'Dark2', 'GnBu',
                    'Greens', 'Greys', 'OrRd', 'Oranges', 'PRGn', 'Paired',
                    'Pastel1', 'Pastel2', 'PiYG', 'PuBu', 'PuBuGn', 'PuOr',
                    'PuRd', 'Purples', 'RdBu', 'RdGy', 'RdPu', 'RdYlBu', 'RdYlGn',
                    'Reds', 'Set1', 'Set2', 'Set3', 'Spectral', 'YlGn', 'YlGnBu',
                    'YlOrBr', 'YlOrRd', 'afmhot', 'autumn', 'binary', 'bone',
                    'brg', 'bwr', 'cool', 'copper', 'flag', 'gist_earth',
                    'gist_gray', 'gist_heat', 'gist_ncar', 'gist_rainbow',
                    'gist_stern', 'gist_yarg', 'gnuplot', 'gnuplot2', 'gray',
                    'hot', 'hsv', 'jet', 'ocean', 'pink', 'prism', 'rainbow',
                    'seismic', 'spectral', 'spring', 'summer', 'terrain', 'winter'] + [n for n,c in Heatmap.my_colormaps]
        for n,c in Heatmap.my_colormaps:
            my_cmap = matplotlib.colors.LinearSegmentedColormap(n,c,256)
            pylab.register_cmap(name=n,cmap=my_cmap)
        arg( '-c','--colormap', type=str, choices = col_maps, default = 'bbcry' )
        arg( '--bottom_c', type=str, default = None,
             help = "Color to use for cells below the minimum value of the scale [default None meaning bottom color of the scale]")
        arg( '--top_c', type=str, default = None,
             help = "Color to use for cells below the maximum value of the scale [default None meaning bottom color of the scale]")
        arg( '--nan_c', type=str, default = None,
             help = "Color to use for nan cells  [default None]")



        """
        arg( '--', type=str, default="average",
             help = "Linkage method for feature clustering [default average]")
        arg( '--slinkage', type=str, default="average",
             help = "Linkage method for sample clustering [default average]")
        """

    def __init__( self, numpy_matrix, sdendrogram, fdendrogram, snames, fnames, fnames_meta, args = None ):
        self.numpy_matrix = numpy_matrix
        self.sdendrogram = sdendrogram
        self.fdendrogram = fdendrogram
        self.snames = snames
        self.fnames = fnames
        self.fnames_meta = fnames_meta
        self.ns,self.nf = self.numpy_matrix.shape
        self.args = args

    def make_legend( self, dmap, titles, out_fn ):
        figlegend = plt.figure(figsize=(1+3*len(titles),2), frameon = False)

        gs = gridspec.GridSpec( 1, len(dmap), wspace = 2.0  )

        for i,(d,title) in enumerate(zip(dmap,titles)):
            legax = plt.subplot(gs[i],frameon = False)
            for k,v in sorted(d.items(),key=lambda x:x[1]):
                rect = Rectangle( [0.0, 0.0], 0.0, 0.0,
                                  facecolor = self.dcols[v%len(self.dcols)],
                                  label = k,
                                  edgecolor='b', lw = 0.0)

                legax.add_patch(rect)
        #remove_splines( legax )
            legax.set_xticks([])
            legax.set_yticks([])
            legax.legend( loc = 2, frameon = False, title = title)
        """
                      ncol = legend_ncol, bbox_to_anchor=(1.01, 3.),
                      borderpad = 0.0, labelspacing = 0.0,
                      handlelength = 0.5, handletextpad = 0.3,
                      borderaxespad = 0.0, columnspacing = 0.3,
                      prop = {'size':fontsize}, frameon = False)
        """
        if out_fn:
            figlegend.savefig(out_fn, bbox_inches='tight')

    def draw( self ):

        rat = float(self.ns)/self.nf
        rat *= self.args.cell_aspect_ratio
        x,y = (self.args.image_size,rat*self.args.image_size) if rat < 1 else (self.args.image_size/rat,self.args.image_size)
        fig = plt.figure( figsize=(x,y), facecolor = 'w'  )

        cm = pylab.get_cmap(self.args.colormap)
        # cm = plt.get_cmap(self.args.colormap)
        bottom_col = [  cm._segmentdata['red'][0][1],
                        cm._segmentdata['green'][0][1],
                        cm._segmentdata['blue'][0][1]   ]
        if self.args.bottom_c:
            bottom_col = self.args.bottom_c
        cm.set_under( bottom_col )
        top_col = [  cm._segmentdata['red'][-1][1],
                     cm._segmentdata['green'][-1][1],
                     cm._segmentdata['blue'][-1][1]   ]
        if self.args.top_c:
            top_col = self.args.top_c
        cm.set_over( top_col )

        if self.args.nan_c:
            cm.set_bad( self.args.nan_c  )

        def make_ticklabels_invisible(ax):
            for tl in ax.get_xticklabels() + ax.get_yticklabels():
                 tl.set_visible(False)
            ax.set_xticks([])
            ax.set_yticks([])

        def remove_splines( ax ):
            for v in ['right','left','top','bottom']:
                ax.spines[v].set_color('none')

        def shrink_labels( labels, n ):
            shrink = lambda x: x[:n/2]+" [...] "+x[-n/2:]
            return [(shrink(str(l)) if len(str(l)) > n else l) for l in labels]


        #gs = gridspec.GridSpec( 4, 2,
        #                        width_ratios=[1.0-fr_ns,fr_ns],
        #                        height_ratios=[.03,0.03,1.0-fr_nf,fr_nf],
        #                        wspace = 0.0, hspace = 0.0 )

        fr_ns = float(self.ns)/max([self.ns,self.nf])
        fr_nf = float(self.nf)/max([self.ns,self.nf])

        buf_space = 0.05
        minv = min( [buf_space*8, 8*rat*buf_space] )
        if minv < 0.05:
            buf_space /= minv/0.05
        metadata_height = self.args.metadata_height if type(snames[0]) is tuple and len(snames[0]) > 1 else 0.000001
        gs = gridspec.GridSpec( 6, 4,
                                width_ratios=[ buf_space, buf_space*2, .08*self.args.fdend_width,0.9],
                                height_ratios=[ buf_space, buf_space*2, .08*self.args.sdend_height, metadata_height, self.args.metadata_separation, 0.9],
                                wspace = 0.0, hspace = 0.0 )

        ax_hm = plt.subplot(gs[23], facecolor = bottom_col  )
        ax_metadata = plt.subplot(gs[15], facecolor = bottom_col  )
        ax_hm_y2 = ax_hm.twinx()

        norm_f = matplotlib.colors.Normalize
        if self.args.log_scale:
            norm_f = matplotlib.colors.LogNorm
        elif self.args.sqrt_scale:
            norm_f = SqrtNorm
        minv, maxv = 0.0, None

        maps, values, ndv = [], [], 0
        if type(snames[0]) is tuple and len(snames[0]) > 1:
            metadata = zip(*[list(s[1:]) for s in snames])
            for m in metadata:
                mmap = dict([(v[1],ndv+v[0]) for v in enumerate(list(set(m)))])
                values.append([mmap[v] for v in m])
                ndv += len(mmap)
                maps.append(mmap)
            dcols = []
            mdmat = np.matrix(values)
            while len(dcols) < ndv:
                dcols += self.dcols
            cmap = matplotlib.colors.ListedColormap(dcols[:ndv])
            bounds = [float(f)-0.5 for f in range(ndv+1)]
            imm = ax_metadata.imshow( mdmat, #origin='lower',
                    interpolation = 'nearest',
                                    aspect='auto',
                                    extent = [0, self.nf, 0, self.ns],
                                    cmap=cmap,
                                    vmin=bounds[0],
                                    vmax=bounds[-1],
                                    )
            remove_splines( ax_metadata )
            ax_metadata_y2 = ax_metadata.twinx()
            ax_metadata_y2.set_ylim(0,len(self.fnames_meta))
            ax_metadata.set_yticks([])
            ax_metadata_y2.set_ylim(0,len(self.fnames_meta))
            ax_metadata_y2.tick_params(length=0)
            ax_metadata_y2.set_yticks(np.arange(len(self.fnames_meta))+0.5)
            ax_metadata_y2.set_yticklabels(self.fnames_meta[::-1], va='center',size=self.args.flabel_size)
        else:
            ax_metadata.set_yticks([])

        ax_metadata.set_xticks([])

        im = ax_hm.imshow( self.numpy_matrix, #origin='lower',
                                interpolation = 'nearest',  aspect='auto',
                                extent = [0, self.nf, 0, self.ns],
                                cmap=cm,
                                vmin=self.args.minv,
                                vmax=self.args.maxv,
                                norm = norm_f( vmin=minv if minv > 0.0 else None, vmax=maxv)
                                )

        #ax_hm.set_ylim([0,800])
        ax_hm.set_xticks(np.arange(len(list(snames)))+0.5)
        if not self.args.no_slabels:
            snames_short = shrink_labels( list([s[0] for s in snames]) if type(snames[0]) is tuple else snames, self.args.max_slabel_len )
            ax_hm.set_xticklabels(snames_short,rotation=90,va='top',ha='center',size=self.args.slabel_size)
        else:
            ax_hm.set_xticklabels([])
        ax_hm_y2.set_ylim([0,self.ns])
        ax_hm_y2.set_yticks(np.arange(len(fnames))+0.5)
        if not self.args.no_flabels:
            fnames_short = shrink_labels( fnames, self.args.max_flabel_len )
            ax_hm_y2.set_yticklabels(fnames_short,va='center',size=self.args.flabel_size)
        else:
            ax_hm_y2.set_yticklabels( [] )
        ax_hm.set_yticks([])
        remove_splines( ax_hm )
        ax_hm.tick_params(length=0)
        ax_hm_y2.tick_params(length=0)
        #ax_hm.set_xlim([0,self.ns])
        ax_cm = plt.subplot(gs[3], facecolor = 'r', frameon = False)
        #fig.colorbar(im, ax_cm, orientation = 'horizontal', spacing = 'proportional', format = ticker.LogFormatterMathtext() )
        cbar = fig.colorbar(im, ax_cm, orientation = 'horizontal', spacing='proportional' if self.args.sqrt_scale else 'uniform' ) # , format = ticker.LogFormatterMathtext() )
        cbar.ax.tick_params(labelsize=args.colorbar_font_size)

        if not self.args.no_sclustering:
            ax_den_top = plt.subplot(gs[11], facecolor = 'r', frameon = False)
            if not self.args.no_plot_sclustering:
                sph._plot_dendrogram( self.sdendrogram['icoord'], self.sdendrogram['dcoord'], self.sdendrogram['ivl'],
                                  self.ns + 1, self.nf + 1, 1, 'top', no_labels=True,
                                  color_list=self.sdendrogram['color_list'] )
            ymax = max([max(a) for a in self.sdendrogram['dcoord']])
            ax_den_top.set_ylim([0,ymax])
            make_ticklabels_invisible( ax_den_top )
        if not self.args.no_fclustering:
            ax_den_right = plt.subplot(gs[22], facecolor = 'b', frameon = False)
            if not self.args.no_plot_fclustering:
                sph._plot_dendrogram(   self.fdendrogram['icoord'], self.fdendrogram['dcoord'], self.fdendrogram['ivl'],
                                    self.ns + 1, self.nf + 1, 1, 'right', no_labels=True,
                                    color_list=self.fdendrogram['color_list'] )
            xmax = max([max(a) for a in self.fdendrogram['dcoord']])
            ax_den_right.set_xlim([xmax,0])
            make_ticklabels_invisible( ax_den_right )

        if self.args.title:
            fig.suptitle(self.args.title,
                         x = 0.5,
                         horizontalalignment = 'center',
                         fontsize = self.args.title_fontsize)
        if not self.args.out:
            plt.show( )
        else:
            fig.savefig( self.args.out, bbox_inches='tight', dpi = self.args.dpi )
            if maps:
                self.make_legend( maps, fnames_meta, self.args.legend_file )



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
        arg( '--legend_file', metavar='LEGEND_FILE', type=str, nargs='?', default=None,
             help= "The output file for the legend of the provided metadata" )

        input_types = [DataMatrix.datatype,DistMatrix.datatype]
        arg( '-t', '--input_type', metavar='INPUT_TYPE', type=str, choices = input_types,
             default='data_matrix',
             help= "The input type can be a data matrix or distance matrix [default data_matrix]" )

        DataMatrix.input_parameters( p )
        DistMatrix.input_parameters( p )
        HClustering.input_parameters( p )
        Heatmap.input_parameters( p )

        self.args  = p.parse_args()

    def check_consistency( self ):
        pass

    def get_args( self ):
        return self.args

if __name__ == '__main__':

    read = ReadCmd( )
    read.check_consistency()
    args = read.get_args()

    if args.input_type == DataMatrix.datatype:
        dm = DataMatrix( args.inp, args )
        if args.out_table:
            dm.save_matrix( args.out_table )

        # print dm.table.axes

        distm = DistMatrix( dm.get_numpy_matrix(), args = args )
        if not args.no_sclustering:
            distm.compute_s_dists()
        if not args.no_fclustering:
            distm.compute_f_dists()
    elif args.input_type == DataMatrix.datatype:
        # distm = read...
        pass
    else:
        pass

    cl = HClustering( distm.get_s_dm(), distm.get_f_dm(), args = args )
    if not args.no_sclustering:
        cl.shcluster()
    if not args.no_fclustering:
        cl.fhcluster()

    hmp = dm.get_numpy_matrix()
    fnames = dm.get_fnames()
    snames = dm.get_snames()
    fnames_meta = snames.names[1:]
    #if not args.no_sclustering or not args.no_fclustering ):

    hmp = cl.get_reordered_matrix( hmp, sclustering = not args.no_sclustering, fclustering = not args.no_fclustering  )
    if not args.no_sclustering:
        snames = cl.get_reordered_sample_labels( snames )
    if not args.no_fclustering:
        fnames = cl.get_reordered_feature_labels( fnames )
    else:
        fnames = fnames[::-1]

    hm = Heatmap( hmp, cl.sdendrogram, cl.fdendrogram, snames, fnames, fnames_meta, args = args )
    hm.draw()
