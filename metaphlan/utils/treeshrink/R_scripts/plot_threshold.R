require("ggplot2")

plotdots = function(d, name=paste("figure",date(),"pdf",sep=".") ) {qplot(1:length(d$V1),d$V1,geom=c("point","line"))+geom_hline(yintercept=0.3170793,color="red",linetype=2)+theme_classic()+xlab("removals")+ylab("delta diameter");ggsave(name,width=4.5,height=3.2);}

args = commandArgs(TRUE)

file = args[1]
outfile = args[2]

d = read.table(file)
plotdots(d,name=outfile)
