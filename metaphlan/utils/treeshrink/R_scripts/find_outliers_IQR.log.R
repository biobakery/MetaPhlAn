args = commandArgs(TRUE)

datafile = args[1]
species = args[2]
genNum = as.numeric(args[3])
k = as.numeric(args[4])

data = read.table(datafile)
#rep = scan("samples.txt")
#d = log(log(c(data[data$V2==species,3],sample(rep,genNum-length(data[data$V2==species,1]),replace=TRUE))))
#d = log(log(data[data$V2==species,3]))
#d = c(d,-d)
d = log(c(data[data$V2==species,3],rep(1,genNum-length(data[data$V2==species,1]))))

q = quantile(d)
IQR = q[4]-q[2]

threshold = exp(q[4]+k*IQR)
data[data$V2==species & data$V3 >= threshold,1]

