
dbmovmf_sample_function <- function(prob){
  sample(1:length(prob), size = 1, replace = FALSE, prob = prob)
}

EMb.E_step <- function(logZt){
  #the following step are necessary to avoid numerical overflows when computing the posterior probabilities
  m = apply(logZt,1,max)
  SumZ = Matrix::rowSums(exp(logZt-m))
  logSumZ =  log(SumZ) + m
  Zt = logZt - logSumZ
  Zt = as(exp(Zt),"dgCMatrix")
  Zt = drop0(Zt)
  Zt
}

SEMb.S_step <- function(logZt){
  n = dim(logZt)[1]
  k = dim(logZt)[2]
  Zt = EMb.E_step(logZt)
  # Perform stochastic column assignement to avoid bad local solutions
  row_partition = apply(Zt,1,dbmovmf_sample_function)
  Z = matrix(0,n,k)
  Z = as(Z,"dgCMatrix")
  Z[cbind(seq_along(row_partition), row_partition)] = 1
  Z
}


CEMb.C_step <- function(logZt){
  n = dim(logZt)[1]
  k = dim(logZt)[2]
  row_partition = apply(logZt,1,which.max)
  Z = matrix(0,n,k)
  Z = as(Z,"dgCMatrix")
  Z[cbind(seq_along(row_partition), row_partition)] = 1
  Z
}

SAEMb.E_step <- function(logZt,iter,max_st_iter){

  if(iter<=max_st_iter)
    Z = SEMb.S_step(logZt)
  else
    Z = EMb.E_step(logZt)
  Z
}


SAEMb.E_step <- function(logZt,iter,max_st_iter){

  if(iter<=max_st_iter)
    Z = SEMb.S_step(logZt)
  else
    Z = EMb.E_step(logZt)
  Z
}

CAEMb.E_step <- function(logZt,iter,max_st_iter){

  if(iter<=max_st_iter)
    Z = SEMb.S_step(logZt)
  else
    Z = CEMb.C_step(logZt)
  Z
}

EMb.update_w <- function(Wt){

  d = dim(Wt)[1]
  k = dim(Wt)[2]
  col_partition = apply(Wt,1,which.max)
  W = matrix(0,d,k)
  W = as(W,"dgCMatrix")
  W[cbind(seq_along(col_partition), col_partition)] = 1
  W
}

SEMb.update_w <- function(Wt){
  d = dim(Wt)[1]
  k = dim(Wt)[2]
  col_partition = apply(Wt,1,dbmovmf_sample_function)
  W = matrix(0,d,k)
  W = as(W,"dgCMatrix")
  W[cbind(seq_along(col_partition), col_partition)] = 1
  W
}

SAEMb.update_w <- function(Wt,iter,max_st_iter){

  if(iter<=max_st_iter)
    W = SEMb.update_w(Wt)
  else
    W = EMb.update_w(Wt)
  W
}




generate_parition <- function(n,k, n_partitions = 1){
    partition = as.integer(sample(as.numeric(1:k),n,replace = TRUE))
    if(n_partitions >1){
		  for(i in 2:n_partitions){
			  partition = rbind(partition,as.integer(sample(as.numeric(1:k),n,replace = TRUE)))
		  }
		  partition = as.matrix(partition)
    }
    partition
  }


nmi <- function(rowcluster,truelabels){
    N_kl = as.matrix(table(rowcluster,truelabels))
    N_k = Matrix::rowSums(N_kl)
    N_l = colSums(N_kl)
    Num = 0
    denum_k = 0
	n = length(rowcluster)
    for(i in 1:nrow(N_kl)){
		denum_k = denum_k + (N_k[i]/n)*log((N_k[i]/n))
		for(j in 1:ncol(N_kl)){
			if(N_kl[i,j]!=0){
				Num = Num +(N_kl[i,j]/n)*log(n*N_kl[i,j]/(N_k[i]*N_l[j]))
			}
		}

    }
    denum_l = 0
    for(j in 1:ncol(N_kl)){
		denum_l = denum_l + (N_l[j]/n)*log((N_l[j]/n))
    }
    resnmi = Num/(sqrt(denum_k*denum_l))
    resnmi
}
