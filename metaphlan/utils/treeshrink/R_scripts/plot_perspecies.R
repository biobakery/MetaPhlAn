args = commandArgs(TRUE)
infile = args[1]
threshold = as.numeric(args[2])
outfile = args[3]

require(ggplot2)

threshold

d = read.table(infile,header=F)
p <- ggplot(d,aes(x=log(V1))) + 
    geom_histogram(bins=100,aes(y=..count..)) +
    #geom_density() + 
    geom_vline(xintercept=log(threshold),color="red",linetype="dashed") + 
    xlab("signature") + theme_classic()
ggsave(outfile,p)
