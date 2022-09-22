threshold <- function(y,e=0.05) {x=y[y>0];qlnorm(p=1-e, sdlog = sd(log(x)), meanlog = mean(log(x)))}

args = commandArgs(TRUE)

datafile = args[1]
e = as.numeric(args[2])

d = read.table(datafile)
threshold(d$V1,e=e)
