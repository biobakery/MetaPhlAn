#threshold <- function(y,e=0.001,epsilon=10^-3) {x=y+epsilon;qlnorm(p=1-e, sdlog = sd(log(x)), meanlog = mean(log(x)))}
threshold <- function(y,e=0.05) {x=y[y>0];qlnorm(p=1-e, sdlog = sd(log(x)), meanlog = mean(log(x)))}

args = commandArgs(TRUE)

datafile = args[1]
e = as.numeric(args[2])

d = read.table(datafile)
t = threshold(d$V1,e=e)
#t = 13.91202
#t = 4.13225
#t = 4.191019
idx=(1:length(d$V1))[d$V1>t]
if (length(idx)){
    idx[length(idx)]
} else {
    0
}

