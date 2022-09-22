args = commandArgs(TRUE)

datafile = args[1]
k = as.numeric(args[2])

d = read.table(datafile)

dl = c(log(d$V1),-log(d$V1))

q = quantile(dl)
IQR = q[4]-q[2]

threshold = exp(q[4]+k*IQR)
threshold
