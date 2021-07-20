# -------------------------------------
# Load packages
# -------------------------------------
uu <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
library(Rcpp)
source("Functions/MTeQTL.R")
source("Functions/HuEst.R")
sourceCpp("Functions/MTeQTL.cpp")
sourceCpp("Functions/HuCpp.cpp")
library(glmnet)
library(QUIC)
library(VIM)
library(laeken)
library(mgcv)
nreps <- 500


# ---------------------------------------------------
# Set model parameters based on slurm array id
# ---------------------------------------------------
model.params1 <- expand.grid(rhoY = rep(c(0.0, 0.1, 0.3, 0.5, 0.7), each = nreps), 
                              numShared_eQTL = 15, 
                              R2 = 0.10)

model.params2 <- expand.grid(rhoY = 0.5, 
                              numShared_eQTL = rep(c(5, 10, 15, 18, 20), each = nreps), 
                              R2 = 0.10)

model.params3 <- expand.grid(rhoY = 0.5, 
                              numShared_eQTL = 15, 
                              R2 = rep(c(0.01, 0.05, 0.10, 0.20, 0.40), each = nreps))
 
model.params <- rbind(model.params1, model.params2, model.params3)
rhoY <- model.params[uu,1]
numShared_eQTL <- model.params[uu,2]
R2 <- model.params[uu,3]
 
# --------------------------------------------------
# Get file paths
# --------------------------------------------------
savename <- paste("rhoY_", 100*rhoY, "_R2_",100*R2 , "_eQTLs_", numShared_eQTL,"_", uu%%nreps + 1,".RDS", sep="")
savefilepath <- paste("Results/", sep="")

# --------------------------------------------------
# Get SNPs from SPATC1L
# --------------------------------------------------
x.raw <- readRDS("genotype_ENSG00000160284.RDS"); 
X <- x.raw

# --------------------------------------------------
# Pruning function based on a correlation cutoff
# --------------------------------------------------
SNP_Prune <- function(X, cor.cutoff) {
  
  # sort by MAF
  keep <- NULL
  Xsorted <- X[,sort(pmin(apply(X, 2, mean), 1 - apply(X, 2, mean)), decreasing=TRUE, index=TRUE)$ix]
  temp <- abs(cor(Xsorted[,1], Xsorted[,-1]))
  k <- dim(Xsorted)[2]
  keepnames <- NULL
  
  while (k > 1) {
    
    keep <- cbind(keep, Xsorted[,1])
    keepnames <- c(keepnames, colnames(Xsorted)[1])
    
    rm <- which(temp >= cor.cutoff)
    Xsorted <- Xsorted[,-c(1, rm+1), drop=FALSE]
    if (dim(Xsorted)[2] == 0) {
      k = 1
    } else {
      temp <- abs(cor(Xsorted[,1], Xsorted[,-1]))
      k <- dim(Xsorted)[2]
      #cat(k, "\n")
    }
  }
  
  return(list("X" = keep, "snpnames" = keepnames))
}

# ----------------------------------------------
# Set dimensions and generate data 
# ----------------------------------------------
q <- 29
ntrain <- 400
nval <- 110
ntest <- 110
n <- 620

set.seed(uu)
train.inds <- sample(1:n, ntrain)
test.inds <- sample(c(1:n)[-c(train.inds)], ntest)
val.inds <- c(1:n)[-c(train.inds, test.inds)]

X.train <- X[train.inds,]
X.val <- X[val.inds,]
X.test <- X[test.inds,]

Xtemp <- SNP_Prune(X, cor.cutoff = .95)
X <- X[,Xtemp$snpnames]
X.train <- X.train[,Xtemp$snpnames]
X.val <- X.val[,Xtemp$snpnames]
X.test <- X.test[,Xtemp$snpnames]
xCor <- cor(X)
p <- dim(X)[2]

# -------------------------------------------
# Functions
# -------------------------------------------
TPR <- function(beta, hatbeta, xCor) {
  q <- dim(beta)[2]
  p <- dim(beta)[1]
  tpcount <- 0
  for (k in 1:q) {
    for (j in 1:p) {
      if (beta[j,k]!=0) {
        if (any(hatbeta[which(abs(xCor[j,]) > .75), k]!=0)) {
          tpcount <- tpcount + 1
        }
      }
    }
  }
  return(tpcount/length(which(beta!=0)))
}

MS <- function(beta, hatbeta) {
  return(sum(hatbeta!=0))
}

named.list = function(...) { 
  l = list(...)
  names(l) = as.character(match.call()[-1])
  l
}


# ------------------------------------
# Generate beta
# ------------------------------------
beta <- matrix(rnorm(p*q, sd = 1), nrow=p)
if (numShared_eQTL > 0) {
  sharedeQTL <- sample(1:p, numShared_eQTL)
  onesVec <- rep(0, p)
  onesVec[sharedeQTL] <- 1
  unShared.mat <- matrix(0, nrow=p, ncol=q)
  for (j in 1:q) {
    unShared.mat[sample(c(1:p)[-sharedeQTL], 20 - numShared_eQTL), j] <- 1
  }
  beta <- beta*(onesVec%*%t(rep(1,q))) + beta*unShared.mat
} else {
  unShared.mat <- matrix(0, nrow=p, ncol=q)
  for (j in 1:q) {
    unShared.mat[sample(c(1:p), 20), j] <- 1
  }
  beta <- beta*unShared.mat
}

# -------------------------------------
# Generate response variables
# -------------------------------------
SigmaYcor <- matrix(0, nrow=q, ncol=q)
SigmaYcor[1:20, 1:20] <- rhoY
SigmaYcor[1:10, 1:10] <- rhoY+.2
diag(SigmaYcor) <- 1

sigmas2 <- rep(0, q)
for (j in 1:q) {
  V1 <- var(rbind(X.train, X.val, X.test)%*%beta[,j])
  sigmas2[j] <- V1*(1 - R2)/(R2)
}

sigmas <- sqrt(sigmas2)
SigmaY <- diag(sigmas)%*%SigmaYcor%*%diag(sigmas)
eoY <- eigen(SigmaY)
SigmaYsqrt <- eoY$vec%*%diag(eoY$val^.5)%*%t(eoY$vec)

Y.train <- X.train%*%beta + matrix(rnorm(q*ntrain), nrow=ntrain)%*%SigmaYsqrt
Y.val <- X.val%*%beta + matrix(rnorm(q*nval), nrow=nval)%*%SigmaYsqrt
Y.test <- X.test%*%beta + matrix(rnorm(q*ntest), nrow=ntest)%*%SigmaYsqrt

SigmaYInv <- solve(SigmaY)
SigmaYInv[which(abs(SigmaYInv) < 1e-8)] <- 0
  
# ---------------------------------------
# Impose missingness at random
# ---------------------------------------
missing.train <- matrix(sample(c(NA, 1), ntrain*q, prob=c(.55, .45), replace=TRUE), nrow=ntrain)
missing.val <- matrix(sample(c(NA, 1), nval*q, prob=c(.55, .45), replace=TRUE), nrow=nval)

Y.GS.train <- Y.train
Y.GS.val <- Y.val
Y.train <- Y.train*missing.train
Y.val <- Y.val*missing.val

#cat("The dimension of the predictor is = ", p, "\n")
alpha <- 2^seq(-1, -9, length=6)
alpha.lasso <- seq(0.1, 1, length=12)
nlambda.Hu <- 100
nlambda <- 30
ntau <- 6
Omega.weight <- matrix(1, nrow=q, ncol=q) 
lambda.lasso <- 10^seq(-2, 2, length=200)
  

# ------------------------------------------------------------
# Prep results for storage 
# ------------------------------------------------------------
R2.lasso.test <- rep(0, q)
R2.lasso.val <- array(0, dim=c(length(alpha.lasso), length(lambda.lasso), q))
R2.lasso.test.1se <- rep(0, q)
R2.inner.test <- array(0, dim=c(length(lambda.lasso), length(alpha.lasso), q))
R2.GS.lasso.test <- rep(0, q)
R2.GS.lasso.val <- array(0, dim=c(length(alpha.lasso), length(lambda.lasso), q))
R2.GS.lasso.test.1se <- rep(0, q)
R2.inner.GS.lasso.val <- array(0, dim=c(length(lambda.lasso), length(alpha.lasso), q))
R2.inner.GS.test <- array(0, dim=c(length(lambda.lasso), length(alpha.lasso), q))
R2.regEM.val <- array(0, dim=c(length(alpha), nlambda, ntau, q))
R2.regEM.test <- array(0, dim=c(length(alpha), nlambda, ntau, q))
TPR.regEM <- array(0, dim=c(length(alpha), nlambda, ntau))
MS.regEM <- array(0, dim=c(length(alpha), nlambda, ntau))
ME.regEM <- array(0, dim=c(length(alpha), nlambda, ntau))
MS.regEM.Omega <- array(0, dim=c(length(alpha), nlambda, ntau))

  
# ----------------------------------------------------------------------
# Tissue-by-tissue elastic net
# ----------------------------------------------------------------------
LassoBeta <- matrix(0, nrow=p, ncol=q)
ptm <- proc.time()

for (respInd in 1:dim(Y.train)[2]) {
  
  Ytemp.train <- as.matrix(Y.train[which(!is.na(Y.train[,respInd])),respInd])
  Xtemp.train <- as.matrix(X.train[which(!is.na(Y.train[,respInd])),])
  R2.inner.val <- array(0, dim=c(length(lambda.lasso), length(alpha.lasso)))
  
  for (alphaInd in 1:length(alpha.lasso)) {
    
    lasso.temp <- glmnet(y = Ytemp.train, x = Xtemp.train, alpha = alpha.lasso[alphaInd], lambda = lambda.lasso)
    lasso.lambda <- lasso.temp$lambda
    
    for (tpInd in 1:length(lasso.temp$lambda)) {
      Y.pred.lasso <- predict(lasso.temp, newx=X.val, s=lasso.temp$lambda[tpInd])
      R2.inner.val[tpInd, alphaInd] <- 1 - sum((Y.val[,respInd] - Y.pred.lasso)^2, na.rm=TRUE)/sum((Y.val[,respInd] - mean(Ytemp.train))^2, na.rm=TRUE)
    }
    
    if (length(lasso.temp$lambda) < length(lambda.lasso)) {
      R2.inner.val[c(length(lasso.temp$lambda)+1):length(lambda.lasso), alphaInd] <- NA
    }
  }
  
  alpha.min.err <- max(which(R2.inner.val >= max(R2.inner.val, na.rm=TRUE) - sd(R2.inner.val)/sqrt(length(R2.inner.val)), arr.ind = TRUE)[,2])
  min.err <- min(which(R2.inner.val[,alpha.min.err] >= max((R2.inner.val - sd(R2.inner.val)/sqrt(length(R2.inner.val)))[,alpha.min.err], na.rm=TRUE)))   
  lasso.temp <- glmnet(y = Ytemp.train, x = Xtemp.train, alpha = alpha.lasso[alpha.min.err], lambda = lambda.lasso)
  Y.pred.lasso <- predict(lasso.temp, newx=X.test, s=lasso.temp$lambda[min.err])
  R2.lasso.test.1se[respInd] <- 1 - sum((Y.test[,respInd] - Y.pred.lasso)^2, na.rm=TRUE)/sum((Y.test[,respInd] - mean(Ytemp.train))^2, na.rm=TRUE)
  
  alpha.min.err <- max(which(R2.inner.val == max(R2.inner.val, na.rm=TRUE), arr.ind = TRUE)[,2])
  min.err <- min(which(R2.inner.val[,alpha.min.err] == max((R2.inner.val)[,alpha.min.err], na.rm=TRUE)))   
  lasso.temp <- glmnet(y = Ytemp.train, x = Xtemp.train, alpha = alpha.lasso[alpha.min.err], lambda = lambda.lasso)
  Y.pred.lasso <- predict(lasso.temp, newx=X.test, s=lasso.temp$lambda[min.err])
  LassoBeta[,respInd] <- coefficients(lasso.temp, s=lasso.temp$lambda[min.err])[-1]
  R2.lasso.test[respInd] <- 1 - sum((Y.test[,respInd] - Y.pred.lasso)^2, na.rm=TRUE)/sum((Y.test[,respInd] - mean(Ytemp.train))^2, na.rm=TRUE)
  cat(respInd, "\n")
  
}

Lasso.time <- proc.time() - ptm
Lasso.R2.test <- R2.lasso.test
Lasso.TPR <- TPR(beta, LassoBeta, xCor)
Lasso.MS <- MS(beta, LassoBeta)

# ----------------------------------------------------------------------
# Gold-standard tissue-by-tissue elastic net
# ----------------------------------------------------------------------  
GSLassoBeta <- matrix(0, nrow=p, ncol=q)
for (respInd in 1:dim(Y.train)[2]) {

  R2.inner.GS.val <- array(0, dim=c(length(lambda.lasso), length(alpha.lasso)))
  
  for (alphaInd in 1:length(alpha.lasso)) {
    
    lasso.temp <- glmnet(y = Y.GS.train[,respInd], x = X.train, alpha = alpha.lasso[alphaInd], lambda = lambda.lasso)
    lasso.lambda <- lasso.temp$lambda
    for (tpInd in 1:length(lasso.temp$lambda)) {
      Y.pred.lasso <- predict(lasso.temp, newx=X.val, s=lasso.temp$lambda[tpInd])
      R2.inner.GS.val[tpInd, alphaInd] <- 1 - sum((Y.val[,respInd] - Y.pred.lasso)^2, na.rm=TRUE)/sum((Y.val[,respInd] - mean(Y.GS.train))^2, na.rm=TRUE)
    }
    
    if (length(lasso.temp$lambda) < length(lambda.lasso)) {
      R2.inner.GS.val[c(length(lasso.temp$lambda)+1):length(lambda.lasso), alphaInd] <- NA
    }
  }
  
  alpha.min.err <- max(which(R2.inner.GS.val >= max(R2.inner.GS.val, na.rm=TRUE) - sd(R2.inner.GS.val)/sqrt(length(R2.inner.GS.val)), arr.ind = TRUE)[,2])
  min.err <- min(which(R2.inner.GS.val[,alpha.min.err] >= max((R2.inner.GS.val - sd(R2.inner.GS.val)/sqrt(length(R2.inner.GS.val)))[,alpha.min.err], na.rm=TRUE)))   
  lasso.temp <- glmnet(y = Y.GS.train[,respInd], x = X.train, alpha = alpha.lasso[alpha.min.err], lambda = lambda.lasso)
  Y.pred.lasso <- predict(lasso.temp, newx=X.test, s=lasso.temp$lambda[min.err])
  R2.GS.lasso.test.1se[respInd] <- 1 - sum((Y.test[,respInd] - Y.pred.lasso)^2, na.rm=TRUE)/sum((Y.test[,respInd] - mean(Y.GS.train[,respInd]))^2, na.rm=TRUE)
  
  alpha.min.err <- max(which(R2.inner.GS.val == max(R2.inner.GS.val, na.rm=TRUE), arr.ind = TRUE)[,2])
  min.err <- min(which(R2.inner.GS.val[,alpha.min.err] == max((R2.inner.GS.val)[,alpha.min.err], na.rm=TRUE)))   
  lasso.temp <- glmnet(y = Y.GS.train[,respInd], x = X.train, alpha = alpha.lasso[alpha.min.err], lambda = lambda.lasso)
  Y.pred.lasso <- predict(lasso.temp, newx=X.test, s=lasso.temp$lambda[min.err])
  GSLassoBeta[,respInd] <- coefficients(lasso.temp, s=lasso.temp$lambda[min.err])[-1]
  R2.GS.lasso.test[respInd] <- 1 - sum((Y.test[,respInd] - Y.pred.lasso)^2, na.rm=TRUE)/sum((Y.test[,respInd] - mean(Y.GS.train[,respInd]))^2, na.rm=TRUE)
  cat(respInd, "\n")
  
}

GS.Lasso.R2.test <- R2.GS.lasso.test
GS.Lasso.TPR <- TPR(beta, GSLassoBeta, xCor)
GS.Lasso.MS <- MS(beta, GSLassoBeta)

# -----------------------------------------------------------
# Hu method ("oracle"/GS for gold-standard)
# -----------------------------------------------------------
R2.Hu.val <- array(0, dim=c(nlambda.Hu, length(alpha), q))
R2.Hu.test <- array(0, dim=c(nlambda.Hu, length(alpha), q))
weight.mat <- matrix(1, nrow=p, ncol=q)
for (alphaInd in 1:length(alpha)) {
  HuFit <- Hu_Estimator(Y = Y.GS.train, X = X.train, alpha = alpha[alphaInd], nlambda = nlambda.Hu, weights = weight.mat, beta.warm = NULL, threshold = .50, max.iters = 1e4, tol = 1e-8, delta.param = .05, Yval = Y.val, Xval = X.val, Ytest = Y.test, Xtest = X.test)
  R2.Hu.val[,alphaInd,] <- HuFit$R2.vec.val
  R2.Hu.test[,alphaInd,] <- HuFit$R2.vec.test
}

inds <- which(apply(R2.Hu.val, c(1,2), mean) == max(apply(R2.Hu.val, c(1,2), mean), na.rm=TRUE), arr.ind = TRUE)
HuFit <- Hu_Estimator(Y = Y.GS.train, X = X.train, alpha = alpha[inds[1,2]], nlambda = nlambda.Hu, weights = weight.mat, beta.warm = NULL, threshold = .50, max.iters = 1e4, tol = 1e-8, delta.param = .05, Yval = Y.val, Xval = X.val, Ytest = Y.test, Xtest = X.test)
GS.HuBeta <- HuFit$beta.list[[inds[1,1]]]
GS.Hu.R2.test <- R2.Hu.test[inds[1,1], inds[1,2], ]
GS.Hu.TPR <- TPR(beta, GS.HuBeta, xCor)
GS.Hu.MS <- MS(beta, GS.HuBeta)

  
# -----------------------------------------------------------
# Hu method (practical version)
# -----------------------------------------------------------
ptm <- proc.time()
R2.Hu.val <- array(0, dim=c(nlambda.Hu, length(alpha), q))
R2.Hu.test <- array(0, dim=c(nlambda.Hu, length(alpha), q))
weight.mat <- sqrt(matrix(max(colSums(!is.na(Y.train)))/colSums(!is.na(Y.train)), nrow=p, ncol=q, byrow=TRUE))
for (alphaInd in 1:length(alpha)) {
  HuFit <- Hu_Estimator(Y = Y.train, X = X.train, alpha = alpha[alphaInd], nlambda = nlambda.Hu, weights = weight.mat, beta.warm = NULL, threshold = .50, max.iters = 1e4, tol = 1e-8, delta.param = .05, Yval = Y.val, Xval = X.val, Ytest = Y.test, Xtest = X.test)
  R2.Hu.val[,alphaInd,] <- HuFit$R2.vec.val
  R2.Hu.test[,alphaInd,] <- HuFit$R2.vec.test
}

inds <- which(apply(R2.Hu.val, c(1,2), mean) == max(apply(R2.Hu.val, c(1,2), mean), na.rm=TRUE), arr.ind = TRUE)
HuFit <- Hu_Estimator(Y = Y.train, X = X.train, alpha = alpha[inds[1,2]], nlambda = nlambda.Hu, weights = weight.mat, beta.warm = NULL, threshold = .50, max.iters = 1e4, tol = 1e-8, delta.param = .05, Yval = Y.val, Xval = X.val, Ytest = Y.test, Xtest = X.test)
Hu.time <- proc.time() - ptm
HuBeta <- HuFit$beta.list[[inds[1,1]]]
Hu.R2.test <- R2.Hu.test[inds[1,1], inds[1,2], ]
Hu.TPR <- TPR(beta, HuBeta, xCor)
Hu.MS <- MS(beta, HuBeta)

# ------------------------------------------------------------
# KNN-imputation
# ------------------------------------------------------------
R2.Hu.val <- array(0, dim=c(nlambda.Hu, length(alpha), q))
R2.Hu.test <- array(0, dim=c(nlambda.Hu, length(alpha), q))
Y.train.k20 <- as.matrix(kNN(data.frame(Y.train),numFun = weightedMean, weightDist = TRUE, k=20)[,1:q])
weight.mat <- matrix(1, nrow=p, ncol=q)
for (alphaInd in 1:length(alpha)) {
  HuFit <- Hu_Estimator(Y = Y.train.k20, X = X.train, alpha = alpha[alphaInd], nlambda = nlambda.Hu, weights = weight.mat, beta.warm = NULL, threshold = .50, max.iters = 1e4, tol = 1e-8, delta.param = .05, Yval = Y.val, Xval = X.val, Ytest = Y.test, Xtest = X.test)
  R2.Hu.val[,alphaInd,] <- HuFit$R2.vec.val
  R2.Hu.test[,alphaInd,] <- HuFit$R2.vec.test
}

inds <- which(apply(R2.Hu.val, c(1,2), mean) == max(apply(R2.Hu.val, c(1,2), mean), na.rm=TRUE), arr.ind = TRUE)
HuFit <- Hu_Estimator(Y = Y.train.k20, X = X.train, alpha = alpha[inds[1,2]], nlambda = nlambda.Hu, weights = weight.mat, beta.warm = NULL, threshold = .50, max.iters = 1e4, tol = 1e-8, delta.param = .05, Yval = Y.val, Xval = X.val, Ytest = Y.test, Xtest = X.test)
HuBeta <- HuFit$beta.list[[inds[1,1]]]
k20.Hu.R2.test <- R2.Hu.test[inds[1,1], inds[1,2], ]
k20.Hu.TPR <- TPR(beta, HuBeta, xCor)
k20.Hu.MS <- MS(beta, HuBeta)

  
# ------------------------------------------------------------
# Penalized-EM method 
# ------------------------------------------------------------
ptm2 <- proc.time()
Y.train.stand <- (Y.train - rep(1, dim(Y.train)[1])%*%t(apply(Y.train, 2, function(x){mean(x, na.rm=TRUE)})))/(rep(1, dim(Y.train)[1])%*%t(apply(Y.train, 2, function(x){sd(x, na.rm=TRUE)})))
weights <- sqrt(matrix(max(colSums(!is.na(Y.train)))/colSums(!is.na(Y.train)), nrow=p, ncol=q, byrow=TRUE))
Omega.warm <- diag(1, q)
beta.warm <- matrix(0, nrow=p, ncol=q)
tpOutput <- MTeQTL_TPGrid(
  Y = Y.train.stand, 
  X = X.train,
  max.iter = 500,
  Omega.warm = Omega.warm,
  beta.warm = beta.warm,
  weights = weights,
  ntau = ntau,
  nlambda = nlambda,
  alpha = alpha,
  Omega.weight = Omega.weight,
  tol = 1e-7)

tau <- tpOutput$tau
lam.mat <- tpOutput$lam.mat
Omega.init <- tpOutput$Omega.init

X.val.stand <- (X.val - tcrossprod(rep(1, dim(X.val)[1]), apply(X.train, 2, mean)))/tcrossprod(rep(1, dim(X.val)[1]), apply(X.train, 2, sd))
X.test.stand <- (X.test - tcrossprod(rep(1, dim(X.test)[1]), apply(X.train, 2, mean)))/tcrossprod(rep(1, dim(X.test)[1]), apply(X.train, 2, sd))
Ytemp.val <- Y.val - tcrossprod(rep(1, dim(Y.val)[1]),colMeans(Y.train, na.rm=TRUE))
Ytemp.test <- Y.test - tcrossprod(rep(1, dim(Y.test)[1]),colMeans(Y.train, na.rm=TRUE))

for (tauInd in 1:ntau) {
  for (alphaInd in 1:length(alpha)) {

    Omega.warm <- Omega.init[[tauInd]]
    beta.warm <- matrix(0, nrow=p, ncol=q)
    Y.stand <- NULL

    for (tpInd in 1:nlambda) {

      temp <- MTeQTL_Glasso(
          Y = as.matrix(Y.train.stand), # matrix of expression (n x K); NAs 
          X = as.matrix(X.train),  # matrix of genotypes (n x p)
          lam1 = lam.mat[tauInd,alphaInd,tpInd], # positive tuning paramter for beta
          lam2 = tau[tauInd], # positive tuning parameter for Omega
          alpha = alpha[alphaInd], 
          max.iter = 500, # maximum number of EM iterations
          tol = 1e-7,
          Omega.warm = Omega.warm, # initializer for Omega
          beta.warm = beta.warm, # initializer for beta
          weights =  weights, 
          Y.stand = Y.stand,
          Omega.weight = Omega.weight)
      
      Omega.warm <- solve(temp$Sigma)
      beta.warm <- temp$beta
      Y.stand <- temp$Y.stand
      
      TPR.regEM[alphaInd, tpInd, tauInd] <- TPR(beta, beta.warm, xCor)
      MS.regEM[alphaInd, tpInd, tauInd] <- MS(beta, beta.warm)
      ME.regEM[alphaInd, tpInd, tauInd] <- sum(diag(tcrossprod(crossprod(beta.warm - beta, crossprod(X.test)), t(beta.warm - beta))))/(110*p*q)
      Omega.warm[which(abs(Omega.warm) < 1e-8)] <- 0
      MS.regEM.Omega[alphaInd, tpInd, tauInd] <- MS(SigmaYInv,  Omega.warm)
      
      # --- get validation set results ------
      Y.pred <-  crossprod(t(X.val.stand), (beta.warm*(tcrossprod(rep(1, p),apply(Y.train, 2, function(x){sd(x,na.rm=TRUE)})))))
      for (hh in 1:q) {
        R2.regEM.val[alphaInd,tpInd,tauInd,hh] <- 1 - sum((Ytemp.val[,hh] - Y.pred[,hh])^2, na.rm=TRUE)/sum(Ytemp.val[,hh]^2, na.rm=TRUE)
      }
      
      # --- get test set results -----------
      Y.pred <-  crossprod(t(X.test.stand), (beta.warm*(tcrossprod(rep(1, p),apply(Y.train, 2, function(x){sd(x,na.rm=TRUE)})))))    
      for (xx in 1:dim(Y.test)[2]) {
        R2.regEM.test[alphaInd,tpInd,tauInd,xx] <-  1 - sum((Ytemp.test[,xx] - Y.pred[,xx])^2, na.rm=TRUE)/sum(Ytemp.test[,xx]^2, na.rm=TRUE)
      }
      
      # ----------------------------------------------------
      # break for-loop if nonsparse and predicting poorly
      # ----------------------------------------------------
      if (tpInd > 5) {
        if ((sum(beta.warm!=0) > .10*p*q) & (tpInd < nlambda) & R2.regEM.val[alphaInd,tpInd,tauInd,] < 
          (max(rowMeans(R2.regEM.val[alphaInd,1:(tpInd-1),tauInd, ]), na.rm=TRUE) - 
            sd(rowMeans(R2.regEM.val[alphaInd,1:(tpInd-1),tauInd, ]))/sqrt(tpInd-1))) {
          R2.regEM.val[alphaInd,(tpInd+1):nlambda,tauInd,] <- NA
          break
        }
      }
      
      if ((sum(Omega.warm > 1e-8)/q^2 > 0.75) & (alphaInd < length(alpha)) & (tauInd < ntau)) {
        R2.regEM.val[(alphaInd+1):length(alpha),(tpInd+1):nlambda,(tauInd+1):ntau,] <- NA
        alphaInd <- length(alpha)
        tauInd <- ntau
        break
      }
    }
  }
}

ours.time <- proc.time() - ptm2

inds <- which(apply(R2.regEM.val, c(1,2,3), mean) == max(apply(R2.regEM.val, c(1,2,3), mean), na.rm=TRUE), arr.ind = TRUE)
ours.TPR <- TPR.regEM[inds[1,1], inds[1,2], inds[1,3]]
ours.MS <- MS.regEM[inds[1,1], inds[1,2], inds[1,3]]
ours.R2.test <- R2.regEM.test[inds[1,1], inds[1,2], inds[1,3],]
ours.MS.Omega <- MS.regEM.Omega[inds[1,1], inds[1,2], inds[1,3]]


saveRDS(named.list(
    Lasso.R2.test,
    Lasso.TPR,
    Lasso.MS,
    GS.Lasso.R2.test,
    GS.Lasso.TPR,
    GS.Lasso.MS,
    GS.Hu.R2.test,
    GS.Hu.TPR,
    GS.Hu.MS,
    Hu.R2.test,
    Hu.TPR,
    Hu.MS,
    k20.Hu.R2.test,
    k20.Hu.TPR,
    k20.Hu.MS,
    ours.TPR,
    ours.MS,
    ours.R2.test,
    ours.MS.Omega,
    ours.time,
    Hu.time,
    Lasso.time
), file = paste(savefilepath, savename, sep=""))
  
  