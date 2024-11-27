#!/usr/bin/env Rscript

### calculate diversity ###

# install required packages

required_pkg <- c("optparse", "ape", "rbiom", "compositions", "BiocManager")
a <- sapply(required_pkg, function(x) {  if (!requireNamespace(x, quietly = TRUE))
  install.packages(x, repos = "http://cran.us.r-project.org")
})
if (! "microbiome" %in% installed.packages()){
  BiocManager::install("microbiome")
}

# accept arguments from command line

library("optparse")

option_list = list(
  
  make_option(c("-f", "--file"), action="store", type="character", default=NULL, 
              help="Merged MetaPhlAn profiles. 
                A table with samples as columns and species as rows is required.",
              metavar="character"),
  
  make_option(c("-o", "--out_directory"), action="store", type="character", default="diversity_analysis",
              help="output directory.
                [default = %default]"),
  
  make_option(c("-p", "--outfile_prefix"), action="store", type="character", default=NULL,
              help="file name prefix of the output distance matrix and log files.
                [default = input file basename]"),
  
  make_option(c("-t", "--tree"), action="store", type="character", default=NULL, 
              help="Full path to the MetaPhlAn species Newick tree.
                Mandatory for computing UniFrac distances."),
  
  make_option(c("-d", "--diversity"), action="store", type="character", default="beta", 
              help="Choose whether to calculate alpha or beta diversity. 
                Options are alpha or beta.
                [default = %default]"),
  
  make_option(c("-m", "--metric"), action="store", type="character", default="bray-curtis", 
              help="Name of the function to use when calculating diversity.
                Options for alpha diversity are richness, shannon, simpson, gini.
                Options for beta diversity are bray-curtis, jaccard, weighted-unifrac, unweighted-unifrac, clr, aitchison.
                [default = %default]"),
  
  make_option(c("-s", "--taxon_separator"), action="store", type="character", default="t__", 
              help="taxon separator used in the input MetaPhlAn table.
                Options are: t__ for MetaPhlAn4 profiles and s__ for MetaPhlAn3 profiles.
                [default = %default]")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop('At least one argument must be supplied (input file).tsv', call.=FALSE)
}

if(! (opt$diversity %in% c('alpha', 'beta'))){
  write(paste0('Method "', opt$diversity, '" not available!'), stdout())
  write(paste0('Available diversity analyses are "alpha" and "beta"'), stdout())
  quit(status = -1)
}

if(opt$diversity =="alpha" & ! (opt$metric %in% c('richness', 'shannon', 'simpson', 'gini'))){
  write(paste0('Method "', opt$metric, '" not available for alpha diversity'), stdout())
  write(paste0('Available alpha-diversity metrics are "richness", shannon", "simpson", "gini".'), stdout())
  quit(status = -1)
}

if(opt$diversity =="beta" & ! (opt$metric %in% c('bray-curtis', 'jaccard', 'weighted-unifrac', 'unweighted-unifrac', 'clr', 'aitchison'))){
  write(paste0('Method "', opt$metric, '" not available for beta diversity'), stdout())
  write(paste0('Available beta-diversity distance functions are "bray-curtis", "jaccard", "weighted-unifrac", "unweighted-unifrac", "clr", "aitchison".'), stdout())
  quit(status = -1)
}

if(! (opt$taxon_separator %in% c('t__', 's__'))){
  write(paste0('Taxon separator "', opt$taxon_separator, '" is not available'), stdout())
  write(paste0('Possible taxon separators are "t__" for MetaPhlAn4 profiles and "s__" for MetaPhlAn3 profiles.'), stdout())
  quit(status = -1)
}

if(is.null(opt$tree) & grepl('unifrac', opt$metric)){
  write(paste0('Selected beta-diversity metric: "', opt$metric, '"'), stdout())
  stop("A  tree is mandatory for computing UniFrac distances. (input tree).nwk", call.=FALSE)
}

for(x in c(opt$file, opt$tree)){
  if(!file.exists(x)){
    stop(paste0('Input file "', x, '" does not exist!'), call.=FALSE)
  }
}

if(is.null(opt$outfile_prefix)){
  outfile_prefix <- basename(opt$file)
  outfile_prefix <- tools::file_path_sans_ext(outfile_prefix)
} else {
  outfile_prefix <- opt$outfile_prefix
}

current_dir <- getwd()
dir.create(file.path(current_dir, opt$out_directory), showWarnings = FALSE)

outfile <- paste(current_dir, opt$out_directory, outfile_prefix, sep="/")

### table preprocessing ###

mpa_table <- read.table(opt$file, comment.char = '#', sep = '\t', header = TRUE, check.names=FALSE)

# check if NCBI id is present

vec <- grepl("ncbi", colnames(mpa_table), ignore.case=TRUE)
if(any(vec)){
  # keep all the columns except the one with NCBI id
  mpa_table <- mpa_table[grep(opt$taxon_separator, mpa_table[,1]), !vec]
} else {
  mpa_table <- mpa_table[grep(opt$taxon_separator, mpa_table[,1]),]
}

if(opt$taxon_separator == "t__"){
  mpa_table[,1] <- gsub(".+\\|t__SGB", "", mpa_table[,1])
} else { 
  mpa_table[,1] <- gsub(".+\\|s__", "", mpa_table[,1])
}

mpa_table[,1] <- gsub("_group$", "", mpa_table[,1])
rownames(mpa_table) <- mpa_table[,1]
mpa_table <- mpa_table[,-1]

# remove samples with all unknowns
removed <- which(colSums(mpa_table) == 0)

if(length(removed)>0){
  
  if(length(removed)==1){
    message = "# WARNING: 1 sample with 100% unknown species was removed from the input table."
  } else {
    message = paste0("# WARNING: ", length(removed), " samples with 100% unknown species were removed from the input table.")
  }
  
  write(message, stdout())
  write(paste(names(removed), collapse='\n'), stdout())
  
  write(message, file=paste0(outfile, '_samples.log'))
  write(paste(names(removed), collapse='\n'), file=paste0(outfile, '_samples.log'), append = TRUE)
  
  # remove samples
  mpa_table <- mpa_table[, -removed]
}

### Data transformation
mpa_table <- mpa_table / 100


### Beta diversity ###

if (opt$diversity == "beta"){
  
  # Bray-Curtis
  if (opt$metric == "bray-curtis"){
    mat <- rbiom::beta.div(as.matrix(mpa_table), method="bray-curtis", weighted=TRUE)
  }
  
  # Jaccard
  if (opt$metric == "jaccard"){
    mat <- rbiom::beta.div(as.matrix(mpa_table), method="jaccard", weighted=FALSE)
  }
  
  # Unifrac
  if (grepl("unifrac", opt$metric)){
    mpa_tree <- ape::read.tree(opt$tree)
    
    if(opt$taxon_separator == "s__"){
      mpa_tree$tip.label <- gsub(".+\\|s__", "", mpa_tree$tip.label)
    }
    
    removed <- setdiff(rownames(mpa_table), mpa_tree$tip.label)
    
    if(length(removed)){
      message = paste0("# WARNING: ", length(removed), " species not present in the tree were removed from the input profile.")
      write(message, stdout())
      write(paste(removed, collapse='\n'), stdout())
      
      write(message, file=paste0(outfile, '_species.log'))
      write(paste(removed, collapse='\n'), file=paste0(outfile, '_species.log'), append = TRUE)
    }
    filt_tree <- ape::keep.tip(mpa_tree, setdiff(rownames(mpa_table), removed))
    filt_mpa_table <- mpa_table[filt_tree$tip.label,]
    
    # check again if after species removal some samples have 0s for all the remaining species, and remove them
    removed <- which(colSums(filt_mpa_table) == 0)
    
    if(length(removed)){
      message = paste0("# WARNING: after removal of species not in the tree, ", length(removed), " samples with 0 abundance of the remaining species were removed from the input table.")
      
      write(message, stdout())
      write(paste(names(removed), collapse='\n'), stdout())
      
      if(file.exists(paste0(outfile, '_samples.log'))){
        write(message, file=paste0(outfile, '_samples.log'), append = TRUE)
      } else { 
        write(message, file=paste0(outfile, '_samples.log'))
      }
      
      write(paste(names(removed), collapse='\n'), file=paste0(outfile, '_samples.log'), append = TRUE)
      
      # remove samples
      filt_mpa_table <- filt_mpa_table[, -removed]
    }
    
    if (opt$metric == "weighted-unifrac"){
      mat <- rbiom::beta.div(as.matrix(filt_mpa_table), tree=filt_tree, method="unifrac", weighted=TRUE)
      
    } else if (opt$metric == "unweighted-unifrac"){
      mat <- rbiom::beta.div(as.matrix(filt_mpa_table), tree=filt_tree, method="unifrac", weighted=FALSE)
    } 
    
  }
  
  # CLR or Aitchison
  if (opt$metric == "clr" || opt$metric == "aitchison"){
    # Centered Log-Ratio
    ait_norm_mpa_table <- compositions::clr(t(mpa_table))
    mat <- t(as.matrix(compositions::as.data.frame.rmult(ait_norm_mpa_table)))

    if (opt$metric == "aitchison"){
      # Aitchison
      mat <- rbiom::beta.div(mat, method="euclidean", weighted=TRUE)
    }
  }
  
  ### Alpha Diversity ###
  
} else if (opt$diversity == "alpha"){
  
  # Richness
  if (opt$metric == "richness"){
    mat <- microbiome::alpha(mpa_table, index = c("richness_observed"))
  }
  
  # Shannon
  if (opt$metric == "shannon"){
    mat <- microbiome::alpha(mpa_table, index = c("diversity_shannon"))
  }
  
  # Simpson
  if (opt$metric == "simpson"){
    mat <- microbiome::alpha(mpa_table, index = c("diversity_gini_simpson"))
  }
  
  # Gini
  if (opt$metric == "gini"){
    mat <- microbiome::alpha(mpa_table, index = c("dominance_gini"))
    row.names(mat) <- colnames(mpa_table)
  }
}

write.table(as.matrix(mat), paste0(outfile, '_', opt$metric, '.tsv'), sep = '\t', quote = FALSE)
