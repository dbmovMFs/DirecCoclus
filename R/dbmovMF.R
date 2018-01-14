#' Diagonal Block Mixture of von Mises-Fisher Dsitributions
#'
#' Provides a simple function to bluid and fit a block mixture of von Mises-Fisher for data co-clustering.
#' @import Matrix
#' @import Bessel
#' @import mclust
#' @param X A sparse or dense data matrix, where rows denote objets lying on a unit-hypersphere and column denote features. Supported format include Matrix, dgCMatrix, dgTMatrix, etc.
#' @param k The number of mixture compoenents (or clusters).
#' @param control A list if additional parameters, see Details.
#' @param ... A list of additional parameters (overriding those specified in control).
#'
#' @details In addition to the parameters described above, the function dbmovMF supports other useful parameters including:
#' * max_iter: the maximum number of iteration
#' * max_st_iter: the maximum number of stochatic iterations. Only useful when performing learning using either SAEMb or CAEMb
#' * n_init: number of time the algorithm will be run with different initializations. The final results correspond to the bust run in terms of log likelihood.
#' * tol: tolerance relatve to the log likelihood to declar convergence.
#' * equal_prop: if true then then all mixture proportion will be consider to be equal, and the mixture propotions parmameters alpha will be ignored.
#' * fit_algo: learning algorithm to use. Supported algorithms include: EMb, CEMb, SEMb, SAEMb and CAEMb.
#' * equal_kappa: if true then all concentration parameters kappa of the differrent mixture components will be equal
#' * kappa_: The value of the concentration parameter. This parameter is only effective when euqal_kappa = TRUE.
#' * row_init: a partition of rows into k cluster to be used for initialization. Could be a vector if n_init = 1 of a matrix of size (n_init x n_rows) if n_init>1. Default = NULL.
#' * col_init: same as row_init but for columns.
#' @md
#' @return returns an object of class dbmovMF including the following attributes:
#' * rowcluster: the partition of rows into k clusters.
#' * colcluster: the partition of columns into k clusters.
#' * kappa_: a vector of concentration parameters.
#' * alpha: a vector of cluster proportions.
#' * alpha: a vector of cluster proportions.
#' * ll: the values of the complete data log likelihood over iterations.
#' * Iter: the total number of iteration needed for convergence.
#'
#' @export
dbmovMF<- function(X,k,control=list(),...){


	#Preparing parameters
	control <- c(control, list(...))

	if(is.null(control$max_iter))
		max_iter = 100
	else
		max_iter = control$max_iter

	if(is.null(control$max_st_iter))
		max_st_iter = round(0.7*max_iter)
	else
		max_st_iter = control$max_st_iter

	if(is.null(control$n_init))
		n_init = 1
	else
		n_init = control$n_init

	if(is.null(control$tol))
		tol = 1e-6
	else
		tol = control$tol

	if(is.null(control$equal_prop))
		equal_prop = FALSE
	else
		equal_prop = control$equal_prop

	if(is.null(control$fit_algo)){
		fit_algo = 'CAEMb'
		print('CAEMb is default fitting algortithm for dbmovMF')
	}
	else
		fit_algo = control$fit_algo

	if(is.null(control$equal_kappa))
		equal_kappa = FALSE
	else{
		equal_kappa = control$equal_kappa
		if(equal_kappa){
			if(is.null(control$kappa_)){
				kappa_ = 100
				print("default value for kappa_ is 100")
			}
			else
				kappa_ = control$kappa_
		}
	}


	if(is.null(control$row_init))
		row_init = NULL
	else{
		row_init = control$row_init
		control$row_init = NULL #to save memory space
	}

	if(is.null(control$col_init))
		col_init = NULL
	else{
		col_init = control$col_init
		control$col_init = NULL #to save memory space
	}


	fit_algorithms <- c("EMb", "CEMb", "SEMb", "SAEMb", "CAEMb")




	#Coherence tests

	if(k>= min(dim(X))){
		stop("more clusters than distinc objects/features (rows/columns)")
	}

	if(equal_kappa){
		if(length(kappa_==1L)){
			if(kappa_<=0)
				stop('kappa_  must be positive')
		}
		else
			stop('kappa_ must be a positive real')
	}

	if(is.na(pmatch(tolower(fit_algo), tolower(fit_algorithms))))
        stop("Invalid fit_algo, please choose one of: SAEMb, CAEMb, EMb, CEMb or SEMb")


	if(n_init>1){
		## Coherence tests for row_init
		if(!is.null(row_init)){
			if(is.matrix(row_init)){
				if(dim(row_init)[1]<n_init)
					stop("Less row partitions than n_init")
				if(dim(row_init)[2]!=nrow(X))
					stop("The length of the row partitions is different from the number of objects")
			}
			else
				stop("Error o_O', row_init is not a matrix")
			}

		## Coherence tests for col_init
		if(!is.null(col_init)){
			if(is.matrix(col_init)){
				if(dim(col_init)[1]<n_init)
					stop("Less column partitions than n_init")
				if(dim(col_init)[2]!=ncol(X))
					stop("The length of the col partitions is different from the number of columns")
			}
			else
				stop("Error o_O', col_init is not a matrix")
		}
	}
	else{
		if(!is.null(row_init)){
			if(length(row_init) != nrow(X))
				stop("The length row_init must be equal to the number of rows of X")
		}
		if(!is.null(col_init)){
			if(length(col_init) != ncol(X))
				stop("The length of the col_init must be equal to the number of columns of X")
		}

	}



	#normalize rows to have unit L2 norm and convert X to column spase matrix
	X = as(X,"dgCMatrix")
	X = X/sqrt(Matrix::rowSums(X*X))

	#useful variables
	n = nrow(X)
	d = ncol(X)
	ll = c(rep(0,max_iter))   #(Expected) classification log-likelihood
	ll_runs = c(rep(0,n_init))
	etp = c(rep(0,max_iter))  # entropy of the latent variable z
	nbIter = max_iter
	vr1 = matrix(c(rep(1,n)),1,n) # 1 vector of dimension (1*n)
	vc1 = matrix(c(rep(1,d)),1,d) # 1 vector of dimension (1*d)
	nu = (d/2)-1                  #useful when dealing with the normalizing constant of the vMF distribution


	#perform one run
	do_one <- function(row_c,col_c){

		#Compute initial binary row-cluster indicator matrix
		Z = matrix(0,n,k)
		Z = as(Z,"dgCMatrix")
		Z[cbind(seq_along(row_c), row_c)] = 1

		#Compute initial binary column-cluster indicator matrix
		W = matrix(0,d,k)
		W = as(W,"dgCMatrix")
		W[cbind(seq_along(col_c), col_c)] = 1



		#row centroids MU^z
		Mu_h = diag(1/sqrt(table(col_c)))
		Mu_h = as(Mu_h,"dgCMatrix")

		alpha = (vr1%*%Z)/n #row cluster proportions, matrix of size (1*k)


		if(!equal_kappa){
			nc = vc1%*%W     #column cluster cardinalities,  matrix of size (1*k)
			#X_w (n*k): compressed X according to column clusters
			X_w = X%*%W
			#X_z_w (k*k): compressed X according to row and column clusters
			X_z_w = Matrix::crossprod(X_w,Z)
			r_bar = Matrix::diag(X_z_w)
			r_bar = r_bar/((alpha[1,]*n)*sqrt(nc[1,]))
			kap = (r_bar*d - r_bar^3)/(1-r_bar^2) #vector containing clusters' kappa parameters
			#kap = rep(10,k)
		}
		else
			kap = rep(kappa_,k)

		KAP = as(Matrix::diag(kap),"dgCMatrix")

		for(iter in 1:max_iter){

			#M-step: starting with M-step is necessary when initializing with skmeans. This also allows us to save computation when evaluating the log-likelihood

			#update W
			Wt = Matrix::crossprod(X,Z)%*%(Mu_h*KAP)

			switch (fit_algo,
                "SAEMb" =  W <- SAEMb.update_w(Wt,iter,max_st_iter),
				"CAEMb" =  W <- SAEMb.update_w(Wt,iter,max_st_iter),  #same as saemb
				"EMb"   =  W <- EMb.update_w(Wt), #same as EMb
				"CEMb"  =  W <- EMb.update_w(Wt), #same as EMb
				"SEMb"  =  W <- SEMb.update_w(Wt),
				           W <- SAEMb.update_w(Wt,iter,max_st_iter) # default CAEMb
			)

			col_partition = apply(W,1,which.max)

			#update Mu_h
			Mu_h = Matrix::diag(1/sqrt(table(col_partition)))
			Mu_h = as(Mu_h,"dgCMatrix")

			#update alpha
			alpha = (vr1%*%Z)/n

			#update kappa
			if(!equal_kappa){
				nc = vc1%*%W     #column cluster cardinalities,  matrix of size (1*k)
				#X_w (n*k): compressed X according to column clusters
				X_w = X%*%W
				#X_z_w (k*k): compressed X according to row and column clusters
				X_z_w = Matrix::crossprod(X_w,Z)
				r_bar = Matrix::diag(X_z_w)
				r_bar = r_bar/((alpha[1,]*n)*sqrt(nc[1,]))
				kap = (r_bar*d - r_bar^3)/(1-r_bar^2) #vector containing clusters' kappa parameters
				KAP = as(Matrix::diag(kap),"dgCMatrix")
			}


			#E-step

			tmpz = rep(0,k)

			logZt = (X%*%W)%*%(Mu_h*KAP)
			if(!equal_kappa)
				tmpz = nu*log(kap) - (nu+1)*log(2*pi) - besselI.nuAsym(kap, nu, k.max = 4, log=TRUE)
			if(!equal_prop)
				tmpz = tmpz + log(alpha[1,])

			logZt = sweep(logZt,2,FUN = "+", tmpz)


			switch (fit_algo,
                "SAEMb" =  Z <- SAEMb.E_step(logZt,iter,max_st_iter),
				"CAEMb" =  Z <- CAEMb.E_step(logZt,iter,max_st_iter),
				"EMb"   =  Z <- EMb.E_step(logZt),
				"CEMb"  =  Z <- CEMb.C_step(logZt),
				"SEMb"  =  Z <- SEMb.S_step(logZt),
				           Z <- CAEMb.E_step(logZt,iter,max_st_iter) # default CAEMb
			)


			#evaluate log-likelihood
			ll[iter] = sum(Z*logZt)
			etp[iter] = -sum(Z@x*log(Z@x))

			if(iter > 1){
				if(abs(ll[iter] - ll[iter-1]) <tol){
					nbIter = iter
					break
				}
			}
		}
		row_partition = apply(Z,1,which.max)
		structure(list(rowcluster = row_partition,colcluster = col_partition,kappa_ = kap,alpha=alpha, ll = ll,iter = nbIter))

	#End of do_one
	}



	#Preparing initial row partition
	if(is.null(row_init)){
		row_c = as.integer(sample(as.numeric(1:k),n,replace = TRUE))
	}
	else{
		if(is.matrix(row_init)){
			row_c = as.integer(row_init[1,])
		}
		else{
			row_c = row_init
		}
	}

	#Preparing initial column partition
	if(is.null(col_init)){
		col_c = as.integer(sample(as.numeric(1:k),d,replace = TRUE))
	}
	else{
		if(is.matrix(col_init)){
			col_c = as.integer(col_init[1,])
		}
		else{
			col_c = col_init
		}
	}


	#Perform run(s)
	run <- do_one(row_c,col_c)
	best <- run$ll[run$iter]
	index_best = 1;

	ll_runs[1] = run$ll[run$iter]
	#in case of multiple runs
	if(n_init >= 2){
		#
		for(i in 2:n_init) {

			#row partition preparation
			if(is.null(row_init)){
				row_c = as.integer(sample(as.numeric(1:k),n,replace = TRUE))
			}
			else{
				row_c = as.integer(row_init[i,])
			}

			#column partition preparation
			if(is.null(col_init)){
				col_c = as.integer(sample(as.numeric(1:k),d,replace = TRUE))
			}
			else{
				col_c = as.integer(col_init[i,])
			}

			rrun <- do_one(row_c,col_c)
			ll_runs[i] = rrun$ll[rrun$iter]
			#update best run is necessary
			if(!is.na(rrun$ll[rrun$iter])){
				if(rrun$ll[rrun$iter] > best){
					run <- rrun
					best <- run$ll[run$iter]
					index_best = i
				}
			}
		}
	#end of runs
	}

	run$index_best = index_best
	run$ll_runs = ll_runs
	run

#end of dbmovMF
}





