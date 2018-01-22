# Directional Co-clustering (an R package)

**Please install the following packages first:**
- devtools
- Bessel 
- mclut
- Matrix

**Install DirecCoclus**
- ```R install_github("DirecCoclus", "dbmovMFs")```

**Usage Example: learning dbmovMFs from the CSTR dataset**
```R
#The CSTR dataset is avaible in folder "Data" of the DirecCoclus repository

#Load CSTR, R.matlab package needed
cstr <- readMat("your local path to cstr.mat")

#Fit dbmovMF to CSTR using SAEMb, see documentation for more details on parameter specification 
resSAEMb <- dbmovMF(cstr$fea,k=4,max_iter = 150,n_init = 10,fit_algo="SAEMb")

#Print the confusion table between the true and estimated clustering
table(resSAEMb$rowcluster,cstr$gnd)

#ARI of row clustering (mclust package needed for this step)
adjustedRandIndex(resSAEMb$rowcluster,cstr$gnd)

#Plot the classification log-likelihood
plot(resSAEMb$ll[1:resSAEMb$iter],type = "l")
```
