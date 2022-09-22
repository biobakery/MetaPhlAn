source("Bootlier/bootlier.R")

args = commandArgs(TRUE)

datafile = args[1]
species = args[2]
genNum = as.numeric(args[3])

data = read.table(datafile)
d = c(data[data$V2==species,3],rep(1,genNum-length(data[data$V2==species,1])))

result = bootstrap.identify.outliers(d)
outliers=result$data.outlier.set
if (length(outliers) > 0){
	s = sort(d)
	q = s[round(length(s)*0.95)]
	t = min(outliers[outliers>q])
	data[data$V2==species & data$V3 >= t,1]
}
