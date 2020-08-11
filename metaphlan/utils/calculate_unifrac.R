#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0)
{
  write('Usage:
          RScript calculate_unifrac.R <merged MetaPhlAn profiles> <Newick tree> <output file> <UniFrac type> <normalization>
        
  <merged MetaPhlAn profiles>\tFull path to the merged MetaPhlAn profiles. A table with samples by columns and species by rows is required.
  <Newick tree>\t\t\tFull path to the MetaPhlAn 3 species Newick tree
  <output file>\t\t\tBasename for output UniFrac distance matrix
  <UniFrac type>\t\tweighted for Weighted UniFrac, unweighted for Unweighted UniFrac (default weighted)
  <normalization>\t\t\tFunction to apply to the data for normalization. Can be none, sinh_sqrt, log10 (default none)
        ', stdout())
  quit(save = 'no')  
}

mpa_infile <- args[1]
tree_file <- args[2]
outfile <- args[3]
weighted <- args[4]=='weighted' || is.na(args[4])
normalize <- args[5] %in% c('sinh_sqrt','log10')


for(x in c(mpa_infile,tree_file)){
if(!file.exists(x)){
  write(paste0('Input file "', x, '" does not exists! Exiting.'), stdout())
  quit(status = -1)
}
}

required_pkg <- c('ape','vegan','rbiom')
a <- sapply(required_pkg, function(x) {  if (!requireNamespace(x, quietly = TRUE))
  install.packages(x)
})


suppressPackageStartupMessages(library(rbiom))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(ape))

mpa_table <- read.table(mpa_infile, comment.char = '#', sep = '\t', header = TRUE)
mpa_table <- mpa_table[grep('s__',mpa_table[,1]),-2]
mpa_table[,1] <- gsub(".+\\|s__", "", mpa_table[,1])
rownames(mpa_table) <- mpa_table[,1]
mpa_table <- mpa_table[,-1]

mpa_tree <- ape::read.tree(tree_file)
mpa_tree$tip.label <- gsub(".+\\|s__", "", mpa_tree$tip.label)

filt_tree <- ape::keep.tip(mpa_tree, intersect(rownames(mpa_table),mpa_tree$tip.label))
filt_mpa_table <- mpa_table[filt_tree$tip.label,] / 100.0

if(normalize == 'log10')
  filt_mpa_table <- log10(1 + filt_mpa_table)
if(normalize == 'sinh_sqrt')
  filt_mpa_table <- asinh(sqrt(filt_mpa_table))
  
rbiom_distmat <- rbiom::unifrac(as.matrix(filt_mpa_table), weighted=weighted, tree=filt_tree)

write.table(as.matrix(rbiom_distmat), outfile,sep = '\t', quote = FALSE)