args = commandArgs(TRUE)
infile = args[1]
outfile = args[2]

require(ggplot2)

d = read.table(infile,header=F)
colnames(d) = c("signature","taxon")
p <- ggplot(d,aes(x=taxon,y=signature)) + geom_col() + 
    scale_y_log10() + 
    xlab("taxon") + ylab("signature") + theme_classic() +
    theme(axis.text.x = element_text(angle=90))
ggsave(outfile,p)
