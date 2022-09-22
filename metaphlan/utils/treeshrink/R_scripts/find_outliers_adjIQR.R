require(robustbase)

args = commandArgs(TRUE)

datafile = args[1]
species = args[2]
genNum = as.numeric(args[3])
k = as.numeric(args[4])

data = read.table(datafile)
s = scan("samples.txt")
d = log(c(data[data$V2==species,3],sample(s,genNum-length(data[data$V2==species,1]),replace=TRUE)))
#d = log(c(data[data$V2==species,3],rep(1,genNum-length(data[data$V2==species,1]))))

box=adjboxStats(d,coef=k)
threshold=exp(box$fence[2])

data[data$V2==species & data$V3 >= threshold,1]

