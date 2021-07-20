MTeQTL_Glasso <- function(

		Y,  # matrix of expression (n x K); NAs 
		X,  # matrix of genotypes (n x p)
		lam1, # positive tuning paramter for beta
		alpha, # tuning parameter in (0,1) 
		lam2, # positive tuning parameter for Omega
		max.iter = 200, # maximum number of EM iterations
		Omega.warm = NULL, # initializer for Omega
		beta.warm = NULL, # initializer for beta
		tol = 1e-7,
		weights = NULL,
		Y.stand = NULL,# set of weights (same dim[2] as X)
    Omega.weight = NULL
    
	){
    
  	# -----------------------------------
  	# prelimiaries
  	# -----------------------------------
  	n <- dim(X)[1]
  	p <- dim(X)[2]
  	q <- dim(Y)[2]

  	# -----------------------------------
  	# get missing indices
  	# -----------------------------------
  	missing.inds <- list(NA)
  	for (k in 1:n) {
  		missing.inds[[k]] <- which(is.na(Y[k,]))
  	}
  	if (max(unlist(lapply(missing.inds, length)))==q) {
  		stop("Please remove row with no observed expression")
  	}

  	# -----------------------------------
  	# Standardize X and Y
  	# -----------------------------------
  	X.stand <- (X - tcrossprod(rep(1, n),apply(X, 2, mean)))/(tcrossprod(rep(1, n), apply(X, 2, sd)))
  	XtX <- crossprod(X.stand)
  	XtXeig <- slanczos(XtX, 1)$values
  	
  	if (is.null(Y.stand)) {
      Y.mean <- colMeans(Y, na.rm=TRUE)
    	Y.stand <- Y - tcrossprod(rep(1, dim(Y)[1]), Y.mean)
    	for (i in 1:n) {
    		if (length(missing.inds[[i]]) > 0) {
    			Y.stand[i,missing.inds[[i]]] <- Y.mean[missing.inds[[i]]]
    		} 
    	}
  	}

  	# -----------------------------------
  	# Get initial values 
  	# -----------------------------------
  	tilde.Sigma <- solve(Omega.warm)
  	Omegatp1 <- Omega.warm
  	tilde.beta <- beta.warm
  	betatp1 <- tilde.beta
  	Delta.up <- tilde.beta
  	Xstandbeta <- crossprod(t(X.stand), betatp1)
    
  	for (kk in 1:max.iter) {

      beta.prev <- betatp1
      Omega.prev <- Omegatp1

  	  # ---------------------------------
  	  # E-step 
  	  # ---------------------------------
  		Et <- Y.stand
  		Vt <- matrix(0, q, q)
  		R <- Et - Xstandbeta
  
      for (k in 1:n) {
        if (length(missing.inds[[k]]) > 0) {
          temp <- solve(tilde.Sigma[-missing.inds[[k]],-missing.inds[[k]], drop=FALSE], tilde.Sigma[-missing.inds[[k]],missing.inds[[k]], drop=FALSE])
          Et[k, missing.inds[[k]]] <- Xstandbeta[k, missing.inds[[k]]] + crossprod(temp, R[k, -missing.inds[[k]]])
          Vt[missing.inds[[k]], missing.inds[[k]]] <- Vt[missing.inds[[k]], missing.inds[[k]]] + tilde.Sigma[missing.inds[[k]], missing.inds[[k]]] - crossprod(tilde.Sigma[-missing.inds[[k]], missing.inds[[k]], drop=FALSE], temp)
        }
      }


  		
  		# ------------------------------------
  		# Conditional M-steps
  		# ------------------------------------
      Omegatp1 <- QUIC(crossprod(Et - Xstandbeta)/n + Vt/n, rho = lam2*Omega.weight, tol = 1e-4, msg = 0, X.init = Omegatp1, W.init = tilde.Sigma)$X		  
		  Sigma <- chol2inv(chol(Omegatp1))
    	betatp1 <- AccProxGrad(Y = Et, X = X.stand, Omega = Omegatp1, 
    						weights = weights, betainit = betatp1, lambda = lam1, alpha =  alpha, 
  						  maxiter = 1e5, tol = 1e-8, XtX = XtX, XtXeig = XtXeig, Qeig = slanczos(Omegatp1, 1)$values)

      Xstandbeta <- crossprod(t(X.stand), betatp1)

      # -----------------------------------
      # Convergence diagnostics
      # ----------------------------------- 
      r1 <- max((betatp1 - tilde.beta)^2)
      r2 <- max((Sigma - tilde.Sigma)^2)
      if(r1 < tol && r2 < tol){
        break
      }

      tilde.Sigma <- Sigma
  	  tilde.beta <- betatp1

  	}
    
  return(list("beta" = betatp1, "Sigma" = Sigma, "Y.stand" = Y.stand 
    #,"obj.func" = obj.func
    ))

} 



MTeQTL_GetGrid2 <- function(
  
  Y,  # matrix of expression (n x K); NAs 
  X,  # matrix of genotypes (n x p)
  lam2, # positive tuning parameter for Omega
  alpha,
  nlambda,
  max.iter = 100, # maximum number of EM iterations
  Omega.warm = NULL, # initializer for Omega
  beta.warm = NULL, # initializer for beta
  weights = NULL, # set of weights (same dim[2] as X)
  Omega.weight = NULL,
  tol = 1e-8
  
){
  
  # --- prelimiaries
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  
  # --- get missing indices
  missing.inds <- list(NA)
  for (k in 1:n) {
    missing.inds[[k]] <- which(is.na(Y[k,]))
  }
  if (max(unlist(lapply(missing.inds, length)))==q) {
    stop("Please remove row with no observed expression")
  }

  # ---- standardize X & Y
  X.stand <- (X - tcrossprod(rep(1, n),apply(X, 2, mean)))/(tcrossprod(rep(1, n), apply(X, 2, sd)))
  XtX <- crossprod(X.stand)
  XtXeig <- slanczos(XtX, 1)$values
  Y.mean <- colMeans(Y, na.rm=TRUE)
  Y.stand <- (Y - rep(1, dim(Y)[1])%*%t(Y.mean))
  for (i in 1:n) {
    if (length(missing.inds[[i]]) > 0) {
      Y.stand[i,missing.inds[[i]]] <- Y.mean[missing.inds[[i]]]
    }
  }
  
  # --- initial values 
  tilde.Sigma <- solve(Omega.warm)
  Omegatp1 <- Omega.warm
  tilde.beta <- beta.warm
  betatp1 <- tilde.beta
  Delta.up <- tilde.beta
  
  for (kk in 1:max.iter) {
    
    # ----- E Step 
    Et <- Y.stand
    Vt <- matrix(0, q, q)
    Xstandbeta <- X.stand%*%tilde.beta
    for (k in 1:n) {
      if (length(missing.inds[[k]]) > 0) {
        temp <- solve(tilde.Sigma[-missing.inds[[k]],-missing.inds[[k]], drop=FALSE], tilde.Sigma[-missing.inds[[k]],missing.inds[[k]], drop=FALSE])
        Et[k,missing.inds[[k]]] <- Xstandbeta[k,missing.inds[[k]]] + crossprod(temp, Et[k,-missing.inds[[k]]] - Xstandbeta[k, -missing.inds[[k]]])
        Vt[missing.inds[[k]], missing.inds[[k]]] <- Vt[missing.inds[[k]], missing.inds[[k]]] + tilde.Sigma[missing.inds[[k]], missing.inds[[k]], drop=FALSE] - crossprod(tilde.Sigma[-missing.inds[[k]], missing.inds[[k]], drop=FALSE], temp)
      }
    }
    
    betatp1 <- tilde.beta
    Omegatp1 <- QUIC(n^(-1)*crossprod(Et  - Xstandbeta) + n^(-1)*Vt,  tol = 1e-4, rho = lam2*Omega.weight, msg = 0, X.init = Omegatp1, W.init = tilde.Sigma)$X
    Sigma <- solve(Omegatp1)
    r1 <- max((betatp1 - tilde.beta)^2)
    r2 <- max((Sigma - tilde.Sigma)^2)
    if (kk > 2) {
      if (r1 < tol  && r2 < tol) {
        break
      }
    }
    tilde.beta <- betatp1
    tilde.Sigma <- Sigma
    
  }

 
  lam.mat <- matrix(Inf, nrow=nlambda, ncol=length(alpha))
  for (alphaInd in 1:length(alpha)) {

    lambda.temp <- 10^seq(2, -4, length=100)
    temp <- -2*crossprod(X.stand, Et%*%Omegatp1)
    check <- matrix(0, nrow=100, ncol=p)
    for (jj in 1:p) {
        for (kk in 1:100) {
          t0 <- pmax(abs(temp[jj,]/n) - lambda.temp[kk]*alpha[alphaInd]*weights[jj,], 0)*sign(temp[jj,]/n)
          check[kk,jj] <- (sqrt(sum(t0^2)) < (1-alpha[alphaInd])*lambda.temp[kk])
          if (check[kk,jj]==0) {
            break
          }
        }
    }

    lam.max <- lambda.temp[max(which(rowSums(check) == p))]
    lam.min <- 0.1*lambda.temp[max(which(rowSums(check) == p))]
    for (kk in 1:nlambda) {
      lam.mat[kk, alphaInd] <- (lam.max^((nlambda - kk)/(nlambda - 1)))*(lam.min^((kk-1)/(nlambda - 1)))
    }

  }

  return(list("lam.mat" = lam.mat, "Sigma"= tilde.Sigma))

}


MTeQTL_GetGrid1 <- function(
  
  Y,  # matrix of expression (n x K); NAs 
  X,  # matrix of genotypes (n x p)
  max.iter = 100, # maximum number of EM iterations
  Omega.warm = NULL, # initializer for Omega
  beta.warm = NULL, # initializer for beta
  weights = NULL, # set of weights (same dim[2] as X)
  Omega.weight = NULL,
  tol = 1e-8

){
  
  # --- prelimiaries
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  
  # --- get missing indices
  missing.inds <- list(NA)
  for (k in 1:n) {
    missing.inds[[k]] <- which(is.na(Y[k,]))
  }

  if (max(unlist(lapply(missing.inds, length)))==q) {
    stop("Please remove row with no observed expression")
  }

  # ---- standardize X & Y
  X.stand <- (X - tcrossprod(rep(1, n),apply(X, 2, mean)))/(tcrossprod(rep(1, n), apply(X, 2, sd)))
  XtX <- crossprod(X.stand)
  XtXeig <- slanczos(XtX, 1)$values
  Y.mean <- colMeans(Y, na.rm=TRUE)
  Y.stand <- (Y - rep(1, dim(Y)[1])%*%t(Y.mean))
  for (i in 1:n) {
    if (length(missing.inds[[i]]) > 0) {
      Y.stand[i,missing.inds[[i]]] <- Y.mean[missing.inds[[i]]]
    }
  }
  
  # --- initial values 
  tilde.Sigma <- solve(Omega.warm)
  Omegatp1 <- Omega.warm
  Xstandbeta <- X.stand%*%beta.warm

  
  for (kk in 1:max.iter) {
    
    # ----- E Step 
    Et <- Y.stand
    Vt <- matrix(0, q, q)
    for (k in 1:n) {
      if (length(missing.inds[[k]]) > 0) {
        temp <- solve(tilde.Sigma[-missing.inds[[k]],-missing.inds[[k]], drop=FALSE], tilde.Sigma[-missing.inds[[k]],missing.inds[[k]], drop=FALSE])
        Et[k,missing.inds[[k]]] <- Xstandbeta[k,missing.inds[[k]]] + crossprod(temp, Et[k,-missing.inds[[k]]] - Xstandbeta[k, -missing.inds[[k]]])
        Vt[missing.inds[[k]], missing.inds[[k]]] <- Vt[missing.inds[[k]], missing.inds[[k]]] + tilde.Sigma[missing.inds[[k]], missing.inds[[k]], drop=FALSE] - crossprod(tilde.Sigma[-missing.inds[[k]], missing.inds[[k]], drop=FALSE], temp)
      }
    }

    # ---- M step 
    Sigma <- n^(-1)*crossprod(Et  - Xstandbeta) + n^(-1)*Vt
    r1 <- 0
    r2 <- max((Sigma - tilde.Sigma)^2)
    if (kk > 2) {
      if (r1 < tol  && r2 < tol) {
        break
      }
    }
    tilde.Sigma <- Sigma

  }
 
  return(list("Sigma" = Sigma))
  
}





MTeQTL_TPGrid <- function(
  Y, 
  X,
  max.iter,
  Omega.warm,
  beta.warm,
  weights,
  ntau,
  nlambda,
  alpha,
  Omega.weight,
  tol) {

  Omega.warm <- diag(1, q)
  beta.warm <- matrix(0, nrow=p, ncol=q)

  temp <- MTeQTL_GetGrid1(
   Y = Y,
   X = X,  # matrix of genotypes (n x p)
   max.iter = 200, # maximum number of EM iterations
   Omega.warm =  Omega.warm, # initializer for Omega
   beta.warm = beta.warm, # initializer for beta
   weights =  weights, # set of weights (same dim[2] as X)
   Omega.weight = Omega.weight,
   tol = tol)

  tau.max <- max(abs(temp$Sigma - diag(diag(temp$Sigma))))
  tau.min <- .1*tau.max
  tau <- rep(0, ntau)
  for (k in 1:ntau) {
    tau[k] <- tau.max^((ntau - k)/(ntau - 1))*tau.min^((k-1)/(ntau - 1))
  }
  lam.mat.ours <- array(0, dim=c(ntau, length(alpha), nlambda))
  Omega.init <- list(NA)

  for (tauInd in 1:ntau) {
   
    temp <- MTeQTL_GetGrid2(
      Y = Y,
      X = X,  # matrix of genotypes (n x p)
      lam2 = tau[tauInd],
      alpha = alpha, 
      nlambda = nlambda,
      max.iter = 200, # maximum number of EM iterations
      Omega.warm =  Omega.warm, # initializer for Omega
      beta.warm = beta.warm, # initializer for beta
      weights =  weights, # set of weights (same dim[2] as X)
      Omega.weight = Omega.weight,
      tol = tol)

   Omega.warm <- solve(temp$Sigma)
   Omega.init[[tauInd]] <- Omega.warm
   lam.mat.ours[tauInd,,] <- t(temp$lam.mat)

  }


  return(list("tau" = tau, "lam.mat" = lam.mat.ours, 
              "Omega.init" = Omega.init))


}




