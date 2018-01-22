# Directional Co-clustering (an R package)

**Please install the following packages first:**
- devtools
- Bessel 
- mclut
- Matrix

**Install DirecCoclus**
- ```install_github("DirecCoclus", "dbmovMFs")```

**Usage Example: learning dbmovMFs from the CSTR dataset**
- Preparing the data
  - Download the CSTR dataset avaible in folder "Data" of DirecCoclus repository
  - Install R.matlab package, necessary for the next step
  - Load data: ```cstr <- readMat("your local path to cstr.mat")```
- Fit dbmovMF to CSTR using SAEMb, see documentation for more details on parameter specification
  - ```res_saemb = dbmovMF(cstr$fea,k=4,max_iter = 150,n_init = 10,fit_algo="SAEMb")```
- ARI of row clustering (mclust package needed for this step): ```adjustedRandIndex(res_saemb$rowcluster,cstr$gnd)```
- Confusion Table: ```table(res_saemb$rowcluster,cstr$gnd)```
- Plot the classification log-likelihood: ```plot(res_saemb$ll[1:res_saemb$iter],type = "l")```



```R
#load CSTR, R.matlab package needed
cstr <- readMat("your local path to cstr.mat")

#Fit dbmovMF to CSTR using SAEMb, see documentation for more details on parameter specification 
res_saemb = dbmovMF(cstr$fea,k=4,max_iter = 150,n_init = 10,fit_algo="SAEMb")```
