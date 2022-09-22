require("ggplot2")

threshold <- function(y,e=0.001,epsilon=10^-3) {x=y+epsilon;qlnorm(p=1-e, sdlog = sd(log(x)), meanlog = mean(log(x)))}

plotdots = function(d, name=paste("figure",date(),"pdf",sep=".") ) {qplot(1:length(d$V1),d$V1,geom=c("point","line"))+geom_hline(yintercept=threshold(d$V1,e=0.01),color="red",linetype=2)+geom_hline(yintercept=threshold(d$V1,e=0.05),color="blue",linetype=2)+theme_classic()+xlab("removals")+ylab("delta diameter");ggsave(name,width=4.5,height=3.2);}

args = commandArgs(TRUE)

file = args[1]
outfile = args[2]

d = read.table(file)
plotdots(d,name=outfile)
