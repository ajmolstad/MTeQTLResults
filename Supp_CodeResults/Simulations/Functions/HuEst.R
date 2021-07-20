# ----------------------------------
# Exact Hu et al estimator 
# ----------------------------------


Hu_Estimator <- function(Y, X, alpha, nlambda = NULL, 
  lambda = NULL, beta.warm = NULL, 
  weights = weights, threshold = .60, max.iters = 1e4, tol = 1e-6, 
  delta.param, Yval, Xval, Ytest, Xtest){
	
	# ---------------------------
	# Preliminaries 
	# ---------------------------
  n <- dim(Y)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]

  Y.m <- apply(Y, 2, function(x){mean(x, na.rm=TRUE)})
  Y.sd <- apply(Y, 2, function(x){sd(x, na.rm=TRUE)})
  X.m <- apply(X, 2, function(x){mean(x, na.rm=TRUE)})
  X.sd <- apply(X, 2, function(x){sd(x, na.rm=TRUE)})

  Y <- (Y - rep(1,n)%*%t(Y.m))/(rep(1,n)%*%t(Y.sd))
  X <- (X - rep(1,n)%*%t(X.m))/(rep(1,n)%*%t(X.sd))
  Yval <- (Yval - rep(1,dim(Yval)[1])%*%t(Y.m))/(rep(1,dim(Yval)[1])%*%t(Y.sd))
  Xval <- (Xval - rep(1,dim(Xval)[1])%*%t(X.m))/(rep(1,dim(Xval)[1])%*%t(X.sd))
  Ytest <- (Ytest - rep(1,dim(Ytest)[1])%*%t(Y.m))/(rep(1,dim(Ytest)[1])%*%t(Y.sd))
  Xtest <- (Xtest - rep(1,dim(Xtest)[1])%*%t(X.m))/(rep(1,dim(Xtest)[1])%*%t(X.sd))

  Y.inner <- Y
  Y.inner[which(is.na(Y))] <- 0
  
  if (is.null(beta.warm)) {
    beta.warm <- matrix(0, nrow=p, ncol=q)
  }
  
  missing.inds <- list(NA)
  for (k in 1:dim(Y)[1]) {
    missing.inds[[k]] <- which(is.na(Y[k,]))
  }
  
  nonmissing.inds <- list(NA)
  for (k in 1:dim(Y)[1]) {
    nonmissing.inds[[k]] <- which(!is.na(Y[k,]))
  }
  
  missing.obs <- list(NA)
  for (k in 1:dim(Y)[2]) {
    missing.obs[[k]] <- which(!is.na(Y[,k]))
  }
  
  missing.indicator <- 1*!is.na(Y)
	tilde.n <- colSums(!is.na(Y))
	weight.mat <- weights
	missing.indicator.tilden <- (rep(1, n)%*%t(tilde.n^(-1)))*missing.indicator
	
  beta.list <- list(NA)


	# ------------------------------------------
	# Compute Lipschitz constant for step size
	# -------------------------------------------
	eig.vec <- rep(0, q)
	for (j in 1:q) {
		eig.vec[j] <- slanczos(tilde.n[j]^(-1)*crossprod(X[missing.obs[[j]],]), 1)$values
	}

	L <- max(eig.vec)
  if (is.null(lambda)) {
  	store <- evalGrad(Yinner = Y.inner, X = X, beta = matrix(0, nrow=p, ncol=q), missingindicator = missing.indicator.tilden)
    temp.seq <- 10^seq(4, -4, length=500)
    check.mat <- matrix(0, nrow=500,ncol=p)
    for (j in 1:p) {
      for (k in 1:500) {
        t0 <-  pmax(abs(store[j,]) - temp.seq[k]*alpha*weight.mat[j,],0)*sign(store[j,])
        check.mat[k,j] <- sqrt(sum(t0^2)) < (1-alpha)*temp.seq[k]
      }
    }

    lambda.max <- temp.seq[max(which(rowSums(check.mat) == p))]
    lam.vec <- rep(0, nlambda)
    lam.min <- (delta.param*lambda.max)
    for (k in 1:nlambda) {
      lam.vec[k] <- lambda.max^((nlambda - k)/(nlambda - 1))*lam.min^((k-1)/(nlambda - 1))
    }
  } else {
    lam.vec <- lambda
    nlambda <- length(lam.vec)
  }

  R2.vec.val <- matrix(0, nrow=length(lam.vec), ncol=q)
  R2.vec.test <- matrix(0, nrow=length(lam.vec), ncol=q)
  R2.vec.both <- matrix(0, nrow=length(lam.vec), ncol=q)

	betakm1 <- beta.warm
	betak <- beta.warm
	betakp1 <- beta.warm
	
	for (kk in 1:length(lam.vec)) {
  	
  	iterating <- TRUE
  	gamma <- 1
  	k.iter <- 1
    obj.old <- evalObj(Yinner = Y.inner, X = X, beta = beta.warm, missingindicator = 	missing.indicator.tilden) + lam.vec[kk]*alpha*sum(weight.mat*abs(beta.warm)) + lam.vec[kk]*(1-alpha)*sum(apply(beta.warm, 1, function(x){sqrt(sum(x^2))}))
    obj.orig <- obj.old
  
  	while (iterating) {
  
  	  theta <- betak + gamma*(betak - betakm1)
  		delta <- theta - L^(-1)*evalGrad(Yinner = Y.inner, X = X, beta = theta, missingindicator = 	missing.indicator.tilden)
  		bar.delta <- pmax(abs(delta) - L^(-1)*lam.vec[kk]*alpha*weight.mat, 0)*sign(delta)
  		if (alpha < 1) {
    		for (j in 1:p) {
    			betakp1[j,] <- max(1 - ((L^(-1)*lam.vec[kk]*(1-alpha))/(sqrt(sum(bar.delta[j,]^2)))), 0)*bar.delta[j,]
    		}
  		} else { 
  		 betakp1 <- bar.delta
  		}
  		
  	  #resid <- max(abs(betakp1 - betak))
      obj <- evalObj(Yinner = Y.inner, X = X, beta = betakp1, missingindicator = 	missing.indicator.tilden) + lam.vec[kk]*alpha*sum(weight.mat*abs(betakp1)) + lam.vec[kk]*(1-alpha)*sum(apply(betakp1, 1, function(x){sqrt(sum(x^2))}))
      if (k.iter %% 1 == 0) {
       # cat(obj, "\n")
      }
      if (k.iter > 10) {
    		if (abs(obj.old - obj) < tol*abs(obj.orig)) {
    			iterating <- FALSE
    		}
      }
  
  		k.iter <- k.iter + 1
  
  		if (k.iter > max.iters) {
  			iterating <- FALSE
  		}
      obj.old <- obj
  		betakm1 <- betak
  		betak <- betakp1
  		gamma <- ((k.iter - 2)/(k.iter + 1))
  	}
    
   cat("# lambda ", kk, ":", sum(betakp1!=0)/length(betakp1), "\n")
    
    beta.list[[kk]] <- tcrossprod(1/X.sd, rep(1,q))*betakp1*tcrossprod(rep(1,p), Y.sd)
    XvalBeta <- Xval%*%betakp1 
    XtestBeta <- Xtest%*%betakp1
    XbothBeta <- rbind(Xval, Xtest)%*%betakp1

    for (xxx in 1:q {

      R2.vec.val[kk,xxx] <- 1 - sum((Yval[,xxx] - XvalBeta[,xxx])^2, na.rm=TRUE)/sum(Yval[,xxx]^2, na.rm=TRUE)
      R2.vec.test[kk,xxx] <- 1 - sum((Ytest[,xxx] - XtestBeta[,xxx])^2, na.rm=TRUE)/sum(Ytest[,xxx]^2, na.rm=TRUE)
      R2.vec.both[kk,xxx] <- 1 - sum((rbind(Yval[,xxx],Ytest[,xxx]) - XbothBeta[,xxx])^2, na.rm=TRUE)/sum(rbind(Yval[,xxx],Ytest[,xxx])^2, na.rm=TRUE)

    }

    if (sum(betakp1!=0)/length(betakp1) > threshold && kk !=length(lam.vec)) {
      for (xxx in 1:q) {
  
        R2.vec.val[(kk+1):nlambda,xxx] <- NA
        R2.vec.test[(kk+1):nlambda,xxx] <- NA
        R2.vec.both[(kk+1):nlambda,xxx] <- NA

      }
      break
    
    }
    weight.val <- colSums(!is.na(Yval))/sum(!is.na(Yval))
    if ((kk > 3) & ((kk + 1) <= nlambda)) {
      if (sum(weight.val*R2.vec.val[kk,]) < .5*max(rowSums((rep(1,length(lam.vec))%*%t(weight.val))*R2.vec.val),na.rm=TRUE)) {
         for (xxx in 1:q) {
            R2.vec.val[(kk+1):nlambda,xxx] <- NA
            R2.vec.test[(kk+1):nlambda,xxx] <- NA
            R2.vec.both[(kk+1):nlambda,xxx] <- NA
        }
      
      break
      }
    }
	}

    named.list <- function(...) { 
      l <- list(...)
      names(l) <- as.character( match.call()[-1] )
      l
    }
    
    
    return(named.list(
                  beta.list,
                  R2.vec.both,
                  R2.vec.val,
                  R2.vec.test))
  
}





Hu_Estimator_GetLambda <- function(Y, X, alpha, nlambda, beta.warm = NULL, weights = weights, threshold = .60, max.iters = 1e4, tol = 1e-6, delta.param){
  
  # ---------------------------
  # Preliminaries 
  # ---------------------------
  n <- dim(Y)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]

  Y.m <- apply(Y, 2, function(x){mean(x, na.rm=TRUE)})
  Y.sd <- apply(Y, 2, function(x){sd(x, na.rm=TRUE)})
  X.m <- apply(X, 2, function(x){mean(x, na.rm=TRUE)})
  X.sd <- apply(X, 2, function(x){sd(x, na.rm=TRUE)})

  Y <- (Y - rep(1,n)%*%t(Y.m))/(rep(1,n)%*%t(Y.sd))
  X <- (X - rep(1,n)%*%t(X.m))/(rep(1,n)%*%t(X.sd))
  Yval <- (Yval - rep(1,dim(Yval)[1])%*%t(Y.m))/(rep(1,dim(Yval)[1])%*%t(Y.sd))
  Xval <- (Xval - rep(1,dim(Xval)[1])%*%t(X.m))/(rep(1,dim(Xval)[1])%*%t(X.sd))
  Ytest <- (Ytest - rep(1,dim(Ytest)[1])%*%t(Y.m))/(rep(1,dim(Ytest)[1])%*%t(Y.sd))
  Xtest <- (Xtest - rep(1,dim(Xtest)[1])%*%t(X.m))/(rep(1,dim(Xtest)[1])%*%t(X.sd))

  Y.inner <- Y
  Y.inner[which(is.na(Y))] <- 0
  
  if (is.null(beta.warm)) {
    beta.warm <- matrix(0, nrow=p, ncol=q)
  }
  
  missing.inds <- list(NA)
  for (k in 1:dim(Y)[1]) {
    missing.inds[[k]] <- which(is.na(Y[k,]))
  }
  
  nonmissing.inds <- list(NA)
  for (k in 1:dim(Y)[1]) {
    nonmissing.inds[[k]] <- which(!is.na(Y[k,]))
  }
  
  missing.obs <- list(NA)
  for (k in 1:dim(Y)[2]) {
    missing.obs[[k]] <- which(!is.na(Y[,k]))
  }
  
  missing.indicator <- 1*!is.na(Y)
  tilde.n <- colSums(!is.na(Y))
  weight.mat <- weights
  missing.indicator.tilden <- (rep(1, n)%*%t(tilde.n^(-1)))*missing.indicator
  
  beta.list <- list(NA)


  # ------------------------------------------
  # Compute Lipschitz constant for step size
  # -------------------------------------------
  eig.vec <- rep(0, q)
  for (j in 1:q) {
    eig.vec[j] <- slanczos(tilde.n[j]^(-1)*crossprod(X[missing.obs[[j]],]), 1)$values
  }

  L <- max(eig.vec)
  store <- evalGrad(Yinner = Y.inner, X = X, beta = matrix(0, nrow=p, ncol=q), missingindicator = missing.indicator.tilden)
  temp.seq <- 10^seq(4, -4, length=500)
  check.mat <- matrix(0, nrow=500,ncol=p)
  for (j in 1:p) {
    for (k in 1:500) {
      t0 <-  pmax(abs(store[j,]) - temp.seq[k]*alpha*weight.mat[j,],0)*sign(store[j,])
      check.mat[k,j] <- sqrt(sum(t0^2)) < (1-alpha)*temp.seq[k]
    }
  }

  lambda.max <- temp.seq[max(which(rowSums(check.mat) == p))]
  lam.vec <- rep(0, nlambda)
  lam.min <- (delta.param*lambda.max)
  for (k in 1:nlambda) {
    lam.vec[k] <- lambda.max^((nlambda - k)/(nlambda - 1))*lam.min^((k-1)/(nlambda - 1))
  }
  
  return(named.list(lam.vec))

}
