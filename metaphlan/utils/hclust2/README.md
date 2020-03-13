Hclust2 is a handy tool for plotting heat-maps with several useful options to produce high quality figures that can be used in publication.

Below is the heatmap produced by Hclust2 on the MetaPhlAn2 abundance profiles of HMP and HMP1-phase2 samples (microbial species and samples are hierarchically clustered).

![Optimized-hmp_profiling_background.png](https://bitbucket.org/repo/gy6GaR/images/1429983578-Optimized-hmp_profiling_background.png)

# Usage #

```
#!python

usage: hclust2.py [-h] [-i [INPUT_FILE]] [-o [OUTPUT_FILE]]
                  [--legend_file [LEGEND_FILE]] [-t INPUT_TYPE] [--sep SEP]
                  [--out_table OUT_TABLE] [--fname_row FNAME_ROW]
                  [--sname_row SNAME_ROW] [--metadata_rows METADATA_ROWS]
                  [--skip_rows SKIP_ROWS] [--sperc SPERC] [--fperc FPERC]
                  [--stop STOP] [--ftop FTOP] [--def_na DEF_NA]
                  [--f_dist_f F_DIST_F] [--s_dist_f S_DIST_F]
                  [--load_dist_matrix_f LOAD_DIST_MATRIX_F]
                  [--load_dist_matrix_s LOAD_DIST_MATRIX_S]
                  [--load_pickled_dist_matrix_f LOAD_PICKLED_DIST_MATRIX_F]
                  [--load_pickled_dist_matrix_s LOAD_PICKLED_DIST_MATRIX_S]
                  [--save_pickled_dist_matrix_f SAVE_PICKLED_DIST_MATRIX_F]
                  [--save_pickled_dist_matrix_s SAVE_PICKLED_DIST_MATRIX_S]
                  [--no_fclustering] [--no_sclustering] [--flinkage FLINKAGE]
                  [--slinkage SLINKAGE] [--dpi DPI] [-l] [--title TITLE] [-s]
                  [--no_slabels] [--minv MINV] [--maxv MAXV] [--no_flabels]
                  [--max_slabel_len MAX_SLABEL_LEN]
                  [--max_flabel_len MAX_FLABEL_LEN]
                  [--flabel_size FLABEL_SIZE] [--slabel_size SLABEL_SIZE]
                  [--fdend_width FDEND_WIDTH] [--sdend_height SDEND_HEIGHT]
                  [--metadata_height METADATA_HEIGHT]
                  [--metadata_separation METADATA_SEPARATION]
                  [--image_size IMAGE_SIZE]
                  [--cell_aspect_ratio CELL_ASPECT_RATIO]
                  [-c {Accent,Blues,BrBG,BuGn,BuPu,Dark2,GnBu,Greens,Greys,OrRd,Oranges,PRGn,Paired,Pastel1,Pastel2,PiYG,PuBu,PuBuGn,PuOr,PuRd,Purples,RdBu,RdGy,RdPu,RdYlBu,RdYlGn,Reds,Set1,Set2,Set3,Spectral,YlGn,YlGnBu,YlOrBr,YlOrRd,afmhot,autumn,binary,bone,brg,bwr,cool,copper,flag,gist_earth,gist_gray,gist_heat,gist_ncar,gist_rainbow,gist_stern,gist_yarg,gnuplot,gnuplot2,gray,hot,hsv,jet,ocean,pink,prism,rainbow,seismic,spectral,spring,summer,terrain,winter,bbcyr,bbcry,bcry}]
                  [--bottom_c BOTTOM_C] [--top_c TOP_C] [--nan_c NAN_C]

optional arguments:
  -h, --help            show this help message and exit
  -i [INPUT_FILE], --inp [INPUT_FILE], --in [INPUT_FILE]
                        The input matrix
  -o [OUTPUT_FILE], --out [OUTPUT_FILE]
                        The output image file [image on screen of not
                        specified]
  --legend_file [LEGEND_FILE]
                        The output file for the legend of the provided
                        metadata
  -t INPUT_TYPE, --input_type INPUT_TYPE
                        The input type can be a data matrix or distance matrix
                        [default data_matrix]

Input data matrix parameters:
  --sep SEP
  --out_table OUT_TABLE
                        Write processed data matrix to file
  --fname_row FNAME_ROW
                        row number containing the names of the features
                        [default 0, specify -1 if no names are present in the
                        matrix
  --sname_row SNAME_ROW
                        column number containing the names of the samples
                        [default 0, specify -1 if no names are present in the
                        matrix
  --metadata_rows METADATA_ROWS
                        Row numbers to use as metadata[default None, meaning
                        no metadata
  --skip_rows SKIP_ROWS
                        Row numbers to skip (0-indexed, comma separated) from
                        the input file[default None, meaning no rows skipped
  --sperc SPERC         Percentile of sample value distribution for sample
                        selection
  --fperc FPERC         Percentile of feature value distribution for sample
                        selection
  --stop STOP           Number of top samples to select (ordering based on
                        percentile specified by --sperc)
  --ftop FTOP           Number of top features to select (ordering based on
                        percentile specified by --fperc)
  --def_na DEF_NA       Set the default value for missing values [default None
                        which means no replacement]

Distance parameters:
  --f_dist_f F_DIST_F   Distance function for features [default correlation]
  --s_dist_f S_DIST_F   Distance function for sample [default euclidean]
  --load_dist_matrix_f LOAD_DIST_MATRIX_F
                        Load the distance matrix to be used for features
                        [default None].
  --load_dist_matrix_s LOAD_DIST_MATRIX_S
                        Load the distance matrix to be used for samples
                        [default None].
  --load_pickled_dist_matrix_f LOAD_PICKLED_DIST_MATRIX_F
                        Load the distance matrix to be used for features as
                        previously saved as pickle file using hclust2 itself
                        [default None].
  --load_pickled_dist_matrix_s LOAD_PICKLED_DIST_MATRIX_S
                        Load the distance matrix to be used for samples as
                        previously saved as pickle file using hclust2 itself
                        [default None].
  --save_pickled_dist_matrix_f SAVE_PICKLED_DIST_MATRIX_F
                        Save the distance matrix for features to file [default
                        None].
  --save_pickled_dist_matrix_s SAVE_PICKLED_DIST_MATRIX_S
                        Save the distance matrix for samples to file [default
                        None].

Clustering parameters:
  --no_fclustering      avoid clustering features
  --no_sclustering      avoid clustering samples
  --flinkage FLINKAGE   Linkage method for feature clustering [default
                        average]
  --slinkage SLINKAGE   Linkage method for sample clustering [default average]


Heatmap options:
  --dpi DPI             Image resolution in dpi [default 150]
  -l, --log_scale       Log scale
  --title TITLE         Title of the plot
  -s, --sqrt_scale      Square root scale
  --no_slabels          Do not show sample labels
  --minv MINV           Minimum value to display in the color map [default
                        None meaning automatic]
  --maxv MAXV           Maximum value to display in the color map [default
                        None meaning automatic]
  --no_flabels          Do not show feature labels
  --max_slabel_len MAX_SLABEL_LEN
                        Max number of chars to report for sample labels
                        [default 15]
  --max_flabel_len MAX_FLABEL_LEN
                        Max number of chars to report for feature labels
                        [default 15]
  --flabel_size FLABEL_SIZE
                        Feature label font size [default 10]
  --slabel_size SLABEL_SIZE
                        Sample label font size [default 10]
  --fdend_width FDEND_WIDTH
                        Width of the feature dendrogram [default 1 meaning
                        100% of default heatmap width]
  --sdend_height SDEND_HEIGHT
                        Height of the sample dendrogram [default 1 meaning
                        100% of default heatmap height]
  --metadata_height METADATA_HEIGHT
                        Height of the metadata panel [default 0.05 meaning 5%
                        of default heatmap height]
  --metadata_separation METADATA_SEPARATION
                        Distance between the metadata and data panels.
                        [default 0.001 meaning 0.1% of default heatmap height]
  --image_size IMAGE_SIZE
                        Size of the largest between width and eight size for

  --cell_aspect_ratio CELL_ASPECT_RATIO
                        Aspect ratio between width and height for the cells of
                        the heatmap [default 1.0]
  -c {Accent,Blues,BrBG,BuGn,BuPu,Dark2,GnBu,Greens,Greys,OrRd,Oranges,PRGn,Paired,Pastel1,Pastel2,PiYG,PuBu,PuBuGn,PuOr,PuRd,Purples,RdBu,RdGy,RdPu,RdYlBu,RdYlGn,Reds,Set1,Set2,Set3,Spectral,YlGn,YlGnBu,YlOrBr,YlOrRd,afmhot,autumn,binary,bone,brg,bwr,cool,copper,flag,gist_earth,gist_gray,gist_heat,gist_ncar,gist_rainbow,gist_stern,gist_yarg,gnuplot,gnuplot2,gray,hot,hsv,jet,ocean,pink,prism,rainbow,seismic,spectral,spring,summer,terrain,winter,bbcyr,bbcry,bcry}, --colormap {Accent,Blues,BrBG,BuGn,BuPu,Dark2,GnBu,Greens,Greys,OrRd,Oranges,PRGn,Paired,Pastel1,Pastel2,PiYG,PuBu,PuBuGn,PuOr,PuRd,Purples,RdBu,RdGy,RdPu,RdYlBu,RdYlGn,Reds,Set1,Set2,Set3,Spectral,YlGn,YlGnBu,YlOrBr,YlOrRd,afmhot,autumn,binary,bone,brg,bwr,cool,copper,flag,gist_earth,gist_gray,gist_heat,gist_ncar,gist_rainbow,gist_stern,gist_yarg,gnuplot,gnuplot2,gray,hot,hsv,jet,ocean,pink,prism,rainbow,seismic,spectral,spring,summer,terrain,winter,bbcyr,bbcry,bcry}
  --bottom_c BOTTOM_C   Color to use for cells below the minimum value of the
                        scale [default None meaning bottom color of the scale]
  --top_c TOP_C         Color to use for cells below the maximum value of the
                        scale [default None meaning bottom color of the scale]
  --nan_c NAN_C         Color to use for nan cells [default None]

```
