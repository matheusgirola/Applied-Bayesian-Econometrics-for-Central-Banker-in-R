################################################################################
# A FAVAR algorithm using Gibbs Sampling, Kalman Filter and Carter Kohn
# following Blake and Mumtaz, Chapter 3
# 
# This is what they called "example4.m"
# 
# Model to estimate: X(it) = b(i)*F(t) + lamda(i)*FFR(t) +  v(it)
#                     Z(t) = c(t) + sum(j = 1 -> P)(B(j)*Z(t-j)) + e(t)
#                     Z(t) = {F(t), FFR(t)}
# X = 40 macroeconomics and financial time series; FFRT = Federal Funds Rate
# from UK, 1970Q1 - 2006Q1
################################################################################

# import packages
library("readxl")
library("pracma")
source("C:/Users/MATHEUS/Desktop/Biblioteca/Trabalhos/R_Trabalhos/Bayesian Methods for Macroeconomists/FunctionsBMM.R")

# import data, names and indexes
names <- read_excel("C:/Users/MATHEUS/Desktop/Biblioteca/Trabalhos/R_Trabalhos/Bayesian Methods for Macroeconomists/Code2017/CHAPTER3/DATA/NAMES.xls",
                    col_names = FALSE)
names <- unlist(names, use.names = FALSE) # transform to character vector
data0 <- read_excel("C:/Users/MATHEUS/Desktop/Biblioteca/Trabalhos/R_Trabalhos/Bayesian Methods for Macroeconomists/Code2017/CHAPTER3/DATA/DATAIN.xls",
                    col_names = FALSE)
data0 <- matrix(unlist(data0), ncol = 40)
colnames(data0) = names
index <- read_excel("C:/Users/MATHEUS/Desktop/Biblioteca/Trabalhos/R_Trabalhos/Bayesian Methods for Macroeconomists/Code2017/CHAPTER3/DATA/INDEX.xls",
                    col_names = FALSE)

# dindex = 1 for series that are log differenced; dindex = 3 differencing without logs
dindex = unlist(index[, 1], use.names = FALSE) 
# index = 1 for 'fast moving' series
index = unlist(index[, 2], use.names = FALSE)

# first difference the data where appropriate
data = matrix(0, nrow = 144, ncol = ncol(data0))
for (i in 1:ncol(data0)){
  if (dindex[i] == 1){
    dat = log(data0[, i])
    dat = diff(dat)*100
  }else if (dindex[i] == 3){
    dat = diff(data0[, i])
  }else{
    dat = data0[2:nrow(data0), i]
  }
  data[ ,i] = dat
}

# Standardise the data
data = standardise(data)

# load policy rate and standardise it
z <- read_excel("C:/Users/MATHEUS/Desktop/Biblioteca/Trabalhos/R_Trabalhos/Bayesian Methods for Macroeconomists/Code2017/CHAPTER3/DATA/BASERATE.xls",
                col_names = "Interest Rate")
z = as.matrix(z)
z = as.matrix(z[2: nrow(z)])
z = standardise(z)

KK = 3          # number of factors
L = 2           # number of lags in the VAR
N = KK + 1      # number of variables in var: K factors plus the interest rate
NN = ncol(data) # size of the panel
T = nrow(data)

#STEP 1: set starting values and priors#########################################
# Get initial guess for the factor via principal components
pmat = extract(data, KK)
pmat = pmat$fac

# state vector S[t-1|t-1]
beta0 = c(pmat[1, ], z[1], matrix(0, nrow = 1, ncol = N))
ns = length(beta0)
P00 = diag(nrow = ns) # P[t-1|t-1]

# Arbitraly set Rn = 1 and Sigma to a indetity matrix to start
rmat = matrix(1, nrow = NN, ncol = 1) # variance of the idiosyncratic component
Sigma = diag(nrow = N)                # Variance of VAR errors

# We use "flat prior" for the factors loadings, variance and VAR. This means we
# wont set a prior distribution, therefore the conditional posterior collapses 
# to a OLS formulae

REPS = 5000; BURN = 4000
irfmat = array(0, dim = c((REPS-BURN), 36, (NN+1))) # store IRF results

# Gibbs Sampling ###############################################################
for (igibbs in 1:REPS){
  print(paste("Gibbs iteration: ", igibbs))
  
  # Step 2: sample factor loadings 
  fload = c()
  floadr = c()
  error = c()
  for (i in 1:NN){
    y = data[, i]
    if (index[i] == 0){
      x = pmat
    }else{
      x = cbind(pmat, z)
    }
    M = solve(t(x)%*%x)%*%t(x)%*%y # H* with flat prior
    V = rmat[i]*solve(t(x)%*%x)  # V* with flat prior
    
    # draw
    ff = M + t(rnorm(ncol(x))%*%cholx(V))
    
    # save
    if (index[i] == 0){
      fload = rbind(fload, t(ff))
      floadr = rbind(floadr, 0)
    }else{
      fload = rbind(fload, t(ff[1:(nrow(ff) - 1)]))
      floadr = rbind(floadr, tail(ff, 1))
    }
    error = cbind(error, (y - (x%*%ff)))
  }
  # for identification top K by K block of fload is identity
  fload[1:KK, 1:KK] = diag(nrow = KK)
  
  # for identification top K by 1 block of floadr is zero
  floadr[24:(24 + KK - 1), 1] = matrix(0, nrow = KK, ncol = 1)
  
  # Step 3: sample variance of the idiosyncratic components from
  # inverse gamma
  rmat = c()
  for (i in 1:NN){
    rmati = IG(0, 0, error[, i])
    rmat = cbind(rmat, rmati)
  }
  
  # Step 4: sample VAR coefficients
  Y = cbind(pmat, z)
  X = cbind(lag0(Y,1), lag0(Y,2), 1)
  Y = Y[2:nrow(Y), ]
  X = X[2:nrow(X), ]
  
  M = matrix(solve(t(X)%*%X)%*%t(X)%*%Y, ncol = 1) # conditional mean
  V = kronecker(Sigma, solve(t(X)%*%X))            # conditional variance
  
  # check stationarity
  chck = -1
  while (chck < 0){
    beta = M + t(rnorm(N*((N*L) + 1))%*%cholx(V)) # VAR coefficients'draws
    S = stability(beta, N, L, 1)
    if (S == 0){
      chck = 10
    }
  }
  beta1 = matrix(beta, nrow = ((N*L) + 1), ncol = N)
  errorsv = Y - X%*%beta1
  
  # Sample VAR covariance
  scale = t(errorsv)%*%errorsv
  Sigma = IWPQ(T, solve(scale))
  
  # Step 5: prepare matrices for the state space
  # Y = H*factors + e
  # factors = MU + F*Factors(-1) + v
  # e ~ N(0, R)
  # v ~ N(0,Q)
  # matrix of factors loadings
  H = matrix(0, nrow = (NN + 1), ncol = ((KK + 1)*L))
  H[1:nrow(fload), 1:(KK + 1)] = cbind(fload, floadr)
  H[(nrow(floadr) + 1), (KK + 1)] = 1
  # matrix R
  R = diag(c(rmat, 0))
  # vector MU
  MU = t(matrix(c(beta1[nrow(beta1), ], rep(0, N*(L-1))), ncol = 1))
  # matrix F
  F = rbind(t(beta1[1:(N*L), ]), diag(1, nrow = N*(L-1), ncol = N*L))
  # matrix Q
  Q = matrix(0, nrow = nrow(F), ncol = nrow(F))
  Q[1:N, 1:N] = Sigma
  
  # Carter and Kohn algorithm to draw the factor
  beta_tt = c()                      # will hold the filtered state variable
  ptt = array(0, dim = c(T, ns, ns)) # will hold its variance
  
  # Step 6a: run Kalman Filter
  i = 1
  x = H # this is no longer data, but a matrix of coefficients
  
  # Prediction
  beta10 = MU + matrix(beta0, nrow = 1)%*%t(F)
  p10 = F%*%P00%*%t(F) + Q
  yhat = t(x%*%t(beta10))
  eta = c(data[i, ], z[i, ]) - yhat
  feta = (x%*%p10%*%t(x)) + R
  
  # Updating
  K = (p10%*%t(x))%*%solve(feta)
  beta11 = t(t(beta10) + K%*%t(eta))
  p11 = p10 - K%*%(x%*%p10)
  beta_tt = rbind(beta_tt, beta11)
  ptt[i, , ] = p11
  
  for (i in 2:T){
    # Prediction
    beta10 = MU + beta11%*%t(F)
    p10 = F%*%p11%*%t(F) + Q
    yhat = t(x%*%t(beta10))
    eta = c(data[i, ], z[i, ]) - yhat
    feta = (x%*%p10%*%t(x)) + R
    
    # Updating
    K = (p10%*%t(x))%*%solve(feta)
    beta11 = t(t(beta10) + K%*%t(eta))
    p11 = p10 - K%*%(x%*%p10)
    beta_tt = rbind(beta_tt, beta11)
    ptt[i, , ] = p11
  }
  # Backward recursion to calculate the mean and variance of the state's distrib
  # vector
  beta2 = matrix(0, nrow = T, ncol = ns) # will hold the state's variable draw
  jv = (1:N); jv1 = (1:KK) # index of state variable to extract
  wa = matrix(rnorm(T*ns), nrow = T, ncol = ns)
  f = F[jv, ]   # F*
  q = Q[jv, jv] # Q*
  mu = MU[jv]   # mu*
  
  # period t
  i = T
  p00 = drop(ptt[i, jv1, jv1]); beta2[i, ] = beta_tt[i, ]
  # draw for beta in period t from N(beta_tt, ptt)
  beta2[i, jv1] = beta_tt[i, jv1] + (wa[i, jv1]%*%cholx(p00))
  # periods t-1 to 1:
  for (i in (T-1):1){
    pt = drop(ptt[i, , ])
    bm = beta_tt[i, ] + t(pt%*%t(f)%*%solve(f%*%pt%*%t(f) + q)%*%
                           t(beta2[(i+1), jv] - mu - beta_tt[i, ]%*%t(f)))
    pm = pt - pt%*%t(f)%*%solve(f%*%pt%*%t(f) + q)%*%f%*%pt
    beta2[i, ] = bm
    beta2[i, jv1] = bm[jv1] + (wa[i, jv1]%*%cholx(pm[jv1, jv1]))
  }
  pmat = beta2[, 1:3] # update factors
  if (igibbs > BURN){
    temp = igibbs - BURN
    
    # Compute IRF
    A0 = cholx(Sigma)
    yhat = matrix(0, nrow = 36, ncol = N)
    vhat = matrix(0, nrow = 36, ncol = N)
    vhat[3, 1:N] = c(0, 0, 0, 1)
    for (i in 3:36){
      yhat[i, ] = c(yhat[(i-1), ], yhat[(i-2), ], 1)%*%
        rbind(beta1[1:(N*L), ], 0) + vhat[i, ]%*%A0
    }
    yhat1 = yhat%*%t(H[ ,1:(KK+1)]) # impulse response for the panel
    irfmat[temp, 1:36, 1:(NN+1)] = yhat1
  }
}

# PLOT IRFS ####################################################################
plotIRF <- function(col_data){
  summary = matrix(0, nrow = 3, ncol = 36)
  rownames(summary) = c("16% percentile", "Median", "68% percentile")
  for (m in 1:36){
    summary[,m] = quantile(irfmat[ , m, col_data], c(0.16, 0.50, 0.84))
  }
  summary = rbind(summary, 0)
  matplot(x = c(1, 36), y = c(min(summary), max(summary)), type =  "n", 
          xlab = "Months after schock", ylab = "Response",
          main = paste("IRF", colnames(data0)[col_data] ))
  matpoints(x = (1:36), y = t(summary), type = "l", lty = 1, lwd = 1.5,)
  legend(x = "bottomright", c("16% percentile", "Median", "84% percentile", "0 line"),
         fill = c("green", "red", "black", "blue"), cex = 0.5)
}

par(mfrow=c(4,5))
for(i in 1:20){
  plotIRF(i)
}

par(mfrow=c(4,5))
for(j in 1:40){
  plotIRF(j)
}
