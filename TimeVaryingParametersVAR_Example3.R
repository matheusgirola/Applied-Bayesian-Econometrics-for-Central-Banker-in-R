################################################################################
# A TVP-VAR algorithm using Gibbs Sampling, Kalman Filter and Carter Kohn
# following Blake and Mumtaz, Chapter 3
# 
# This is what they called "example3.m"
# 
# Model to estimate: Y(t) = Beta(t)*X(t) + e1
#                    Beta(t) = mu + F*Beta(t-1) + e2
# T = dlog(GDP Growth), dlog(Inflation) and Federal Funds Rate
################################################################################

# import packages
library("readxl")
source("C:/Users/MATHEUS/Desktop/Biblioteca/Trabalhos/R_Trabalhos/Bayesian Methods for Macroeconomists/FunctionsBMM.R")

# import data
data0 <- read_excel("C:/Users/MATHEUS/Desktop/Biblioteca/Trabalhos/R_Trabalhos/Bayesian Methods for Macroeconomists/Code2017/CHAPTER3/DATA/usdata.xls",
                    col_names = c("GDP Growth", "Inflation", "Federal Funds Rate"))
data = as.matrix(data0)/100
N = ncol(data)
L = 2 # number of lags in the VAR

# Take lags
Y = data
X = cbind(lag0(Y,1), lag0(Y,2), 1)
Y = Y[3:nrow(Y), ]
X = X[3:nrow(X), ]

# Step 1: set starting values and priors using a pre-sample of 10 years ########
T0 = 40 # is quartely data, so 4*10 years
y0 = Y[1:T0, ]
x0 = X[1:T0, ]
b0 = solve(t(x0)%*%x0)%*%t(x0)%*%y0
e0 = y0 - x0%*%b0
sigma0 = (t(e0)%*%e0)/T0
V0 = kronecker(sigma0, solve(t(x0)%*%x0)) # p[0|0]

# priors for the variance of the transition equation
Q0 = V0*T0*(3.5e-04) # p[0|0]*T0*tau -> prior for the variance of the transition
                     # equation error
P00 = V0            # p[0|0] -> variance of the initial state vector. Variance
                     # of state variable p[t-1|t-1]
beta0 = t(matrix(b0, ncol = 1)) # initial state vector. State variable b[t-1|t-1]

# Step 2: Gibbs Sampling #######################################################
# Initialize
Q = Q0
R = sigma0

# remove initial sample
Y = Y[(T0+1):nrow(Y), ]
X = X[(T0+1):nrow(X), ]
T = nrow(X)

REPS = 5000 # 110000
BURN = 4000 # 109000
out1 = array(0, dim = c((REPS-BURN), T, ncol(beta0)))
out2 = array(0, dim = c((REPS-BURN), N, N))
out3 = array(0, dim = c((REPS-BURN), N*((N*L)+1), N*((N*L)+1)))

for (igibbs in 1:REPS){
  print(paste("Gibbs iteration: ", igibbs))
  
  # Step 2a: Kalman Filter #####################################################
  # Step up matrices
  ns =  ncol(beta0)
  F = diag(nrow = ns)                # fixed
  mu = 0                             # fixed
  beta_tt = c()                      # holds filtered state variable
  ptt = array(0, dim = c(T, ns, ns)) # holds filtered state variable variance
  beta11 = beta0
  p11 = P00
  
  # Run the filter
  for (i in 1:T){
    x =  kronecker(diag(nrow = N), matrix(X[i, ], nrow = 1))
    
    # Prediction
    beta10 = mu + beta11%*%t(F)
    p10 = F%*%p11%*%t(F) + Q
    yhat = t(x%*%t(beta10))
    eta = Y[i, ] - yhat
    feta = (x%*%p10%*%t(x)) + R
    
    # Updating
    K = (p10%*%t(x))%*%solve(feta)
    beta11 = t(t(beta10) + K%*%t(eta))
    p11 = p10 - K%*%(x%*%p10)
    ptt[i, , ] = p11
    beta_tt = rbind(beta_tt, beta11)
  }
  # End of Kalman Filter
  # Step 2b: Backward recursion ################################################
  # Step to compute the mean variance of the state space vector's distribution
  chck = -1
  
  # While loop to ensure VAR stable each point in time
  while (chck < 0){
    beta2 = matrix(0, nrow = T, ncol = ns) # holds the draw of the state variable
    wa = matrix(rnorm(T*ns), nrow = T, ncol = ns)
    error = matrix(0, nrow = T, ncol = N)
    roots = matrix(0, nrow = T, ncol = 1)
    
    # period t
    i = T 
    p00 = drop(ptt[i, , ])
    # draw from beta in period t from N(beta_tt, ptt)
    beta2[i, ] = beta_tt[i, ] + (wa[i, ]%*%chol(p00))
    # var residuals (compute in the same loop for convenience)
    error[i, ] = Y[i, ] - X[i, ]%*%matrix(beta2[i, ], nrow = ((N*L) + 1), ncol = N)
    # check stability at ith time, roots[i] = 1 if stability is violated
    roots[i] = stability(t(beta2[i, ]), N, L, 1)
    
    # period t-1 to 1
    for (i in (T-1):1){
      pt = drop(ptt[i, , ])
      bm = beta_tt[i, ] + t(pt%*%t(F)%*%solve(F%*%pt%*%t(F) + Q)%*%
                              +t(beta2[i + 1,] - mu - beta_tt[i, ]%*%t(F)))
      pm = pt - pt%*%t(F)%*%solve(F%*%pt%*%t(F) + Q)%*%F%*%pt
      beta2[i, ] = bm + (wa[i, ]%*%chol(pm))
      error[i, ] = Y[i, ] - X[i, ]%*%matrix(beta2[i, ], nrow = ((N*L) + 1), ncol  = N)
      roots[i] = stability(t(beta2[i, ]), N, L, 1)
    }
    if (sum(roots) == 0){
      chck = 1
    }
  }
  # End of backward recursion
  # Step 3: sample Q from the IW distribution ##################################
  errorq = diff(beta2)
  scaleQ = (t(errorq)%*%errorq) + Q0
  Q = IWPQ((T + T0), solve(scaleQ))
  
  # Step 4: sample R from the IW distribution ##################################
  scaleR = (t(error)%*%error)
  R = IWPQ(T, solve(scaleR))
  
  if(igibbs > BURN){
    temp = igibbs - BURN
    out1[temp, (1:T), ] = beta2
    out2[temp, (1:N), (1:N)] = R
    out3[temp, (1:(N*((N*L) + 1))), (1:(N*((N*L) + 1)))] = Q
  }
}

# Compute IRF to a policy shock using sign restrictions ########################
horz = 40 # impulse response horizon
irfmat = array(0, dim = c((REPS - BURN), T, horz, N)) # matrix to save IR

for (i in 1:(REPS - BURN)){
  print(paste("IRF Step: ", i))
  sigma = drop(out2[i, , ])
  
  # sign restrictions
  chck = -1
  while (chck < 0){
    K = matrix(rnorm(N*N), nrow = N, ncol = N)
    QQ = getqr(K)
    A0hat = chol(sigma)
    A0hat1 = QQ%*%A0hat # candidate draw
    
    for (m in 1:N){
      # check signs in each row
      e1 = A0hat1[m, 1] < 0 # Response of Y
      e2 = A0hat1[m, 2] < 0 # Response of P
      e3 = A0hat1[m, 3] > 0 # Response of R
      
      if ((e1 + e2 + e3) == 3){
        MP = A0hat1[m, ]
        chck = 10
      }else{
        # check signs but reverse them
        e1 = -A0hat1[m, 1] < 0 # Response of Y
        e2 = -A0hat1[m, 2] < 0 # Response of P
        e3 = -A0hat1[m, 3] > 0 # Response of R
        if ((e1 + e2 + e3) == 3){
          MP = -A0hat1[m, ]
          chck = 10
        }
      }
    }
  }
  # re-shuffle rows of A0hat1 and insert MP in the first row
  A0x = c() # will holds rows of A0hat1 not equal to MP
  for (m in 1:N){
    ee = sum(abs(A0hat1[m, ]) == abs(MP))
    if (ee == 0){
      A0x = rbind(A0x, A0hat1[m, ])
    }
  }
  A0new = rbind(A0x, MP) # A0 matrix to be used in the impulse response
  shock = c(0, 0, 1)
  for (j in 1:T){
    btemp = drop(out1[i, j, ])
    btemp = matrix(btemp, nrow = ((N*L) + 1), ncol = N)
    zz = irfsim(btemp, N, L, A0new, shock, (horz + L))
    zz = zz/matrix(zz[1,3], nrow = horz, ncol = N)
    irfmat[i, j, , ] = zz
  }
}

# Plot IRF time varying plot
library(plot3D)
TT = seq(1964.75, 2010.5, 0.25)
HH = 0:(horz - 1)

irf3D <- function(x, horz){
  irf = matrix(0, nrow = length(TT), ncol = horz)
  for (h in 1:horz){
    # in matlab, median(A, 1) -> returns median of COLUMN, to do the samne in R,
    # in the apply function we input margin =2, which corresponds to column
    irf[, h] = drop(apply(irfmat[ , , h, x], 2, median, na.rm=T))
  }
  return(irf)
}

irf1 = irf3D(1, horz)
irf2 = irf3D(2, horz)
irf3 = irf3D(3, horz)   

persp3D(TT, HH, irf1, phi = 25, theta = 135, cex.axis = 0.5, cex.lab = 0.8,
        xlab = "Year", ylab = "Horizon", zlab = "IRF", ticktype = "detailed") 
persp3D(TT, HH, irf2, phi = 25, theta = 135, cex.axis = 0.5, cex.lab = 0.8,
        xlab = "Year", ylab = "Horizon", zlab = "IRF", ticktype = "detailed")
persp3D(TT, HH, irf3, phi = 25, theta = 135, cex.axis = 0.5, cex.lab = 0.8,
        xlab = "Year", ylab = "Horizon", zlab = "IRF", ticktype = "detailed")
