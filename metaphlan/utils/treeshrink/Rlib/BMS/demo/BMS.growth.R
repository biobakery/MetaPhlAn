
data(datafls)

#now do an MC3 sampling over the growth data 'datfls' with 1000 burn-ins,
#9000 iterations (ex burn-ins), 
#and retaining the best 100 models (besides overall MC3 frequencies)
invisible(readline("hit <Return> to do estimate a short BMA MC3 sampling chain."))
 mfls =bms(X.data=datafls,burn=1000,iter=9000,nmodel=100,user.int=T)
#The console printout shows some information on the resulting bma object 'mfls'.
#The chart details the prior and posterior model size (above), and the likelihoodos and MCMC frequencies of the best 100 models (below).
#Setting user.int=FALSE suppresses his kind of output.\n\n ")
invisible(readline("Hit <Return> for some low-level functions."))

# some results:
 summary(mfls)
# summary.bma() shows basic aggregate results from MC3 sampling 

invisible(readline("Hit <Return> for functions on the coefficients."))
 coef(mfls)
 #based on MCMC frequencies, coef.bma shows 
 # the posterior inclusion probabilities (column 1),
 # the unconditional expecteded value of coefficients (column 2),
 # their standard deviations (column 3),
 # the percentage of times the coefficents had a positive sign 
 # (conditional on inclusion, column 4) 

invisible(readline("Hit <Return> for a different version.")) 
 coef(mfls,exact=T,std.coefs=T) 
 #this is similar to coef.bma(), 
 # however here the numbers are based on the exact marginal 
 # likelihoods of the best (100) models drawn. 
 # Moreover the coefficents are shown in standardized form.

invisible(readline("Hit <Return> for other low-level commands."))  
 pmp.bma(mfls) 
  #post. model probs. for top 100 models based on MCMC freqs and likelihoods

beta.draws.bma(mfls[1:3])
 # show the estimates for the best 3 models 
 # the column names are the inclusion vectors in hex-code 
 # (e.g. 101 for inclusion of variables 1 and 3 is "5" in hexciode ) 

mfls[3]$topmod 
 #show the third-best model


invisible(readline("Hit <Return> for some plots ")) 
invisible(par(ask=TRUE))
  density(mfls,reg="BlMktPm")
#plot density for regressor "BlMktPm"

 image(mfls[1:20],FALSE) 
 #plot signs (pos=blue, neg=red, not included=white) for best 20 models


 plotModelsize(mfls,exact=TRUE) 
 #plot prior and posterioro model size based on exact likelhoods 
 # of best (100) models



 plot(mfls)
 #a combined plot
 
invisible(par(ask=FALSE))