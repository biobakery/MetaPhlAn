require("fdrtool")

threshold <- function(y,e=0.05) {x=c(log(y),-log(y)); exp(qhalfnorm(p=1-e, theta=sd2theta(sd(x))))}

args = commandArgs(TRUE)

datafile = args[1]
e = as.numeric(args[2])

d = read.table(datafile)
threshold(d$V1,e=e)
