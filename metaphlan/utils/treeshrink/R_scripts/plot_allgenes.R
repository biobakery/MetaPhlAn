args = commandArgs(TRUE)
infile = args[1]
outfile = args[2]

require(ggplot2)

d = read.table(infile,header=F)
colnames(d) = c("signature")
p <- ggplot(d,aes(x=log(signature))) + 
    geom_histogram(bins=100) + 
    xlab("signature") + ylab("density") + theme_classic() +
    theme(axis.text.x = element_text(angle=90))
ggsave(outfile,p)
