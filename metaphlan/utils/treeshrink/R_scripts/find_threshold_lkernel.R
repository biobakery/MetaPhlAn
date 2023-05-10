threshold <- function(y,e=0.05) {x=y[y>0];exp(quantile.density(density(log(x),adjust=1),p=1-e));}


args = commandArgs(TRUE)

#libpath = paste(args[1],"/Rlib",sep="")
libpath = file.path(args[1],"Rlib")
datafile = args[2]
e = as.numeric(args[3])

suppressMessages(require(BMS,lib.loc=libpath))

d = read.table(datafile)
threshold(d$V1,e=e)
