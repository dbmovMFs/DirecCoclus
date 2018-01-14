# Directional Co-clustering (an R package)

###### Please install the folloing packages first:
-devtools
-Bessel 
-mclut
-Matrix

###### Install DirecCoclus
install_github("DirecCoclus", "dbmovMFs")

###### Download the cstr dataset avaible in the github repository DirecCoclus (In folder Data)
###### Load data, you need the R.matlab package
###### ctsr data, 4 clusters
cstr <- readMat("your local path to cstr.mat")

###### fit dbmovMF to the data, see documentation for parameter specification
res_saemb = dbmovMF(cstr$fea,k=4,max_iter = 150,n_init = 10,fit_algo="SAEMb")

###### you need mclust package for ARI
adjustedRandIndex(res_saemb$rowcluster-1,cstr$gnd)

###### plot le classification log-likelihood
plot(res_saemb$ll[1:res_saemb$iter],type = "l")
