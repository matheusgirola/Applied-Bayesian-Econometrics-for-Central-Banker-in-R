################################################################################
# VAR estimation via Gibbs sampling with conditional forecasting
# following Blake and Mumtaz, Chapter 2
# 
# This is what they called "example8.m"
# 
# Model to estimate: VAR(2) Y(t) = alpha + B(1)Y(t-1) + B(2)Y(t-2) + v(t)
# Y = quartely annual US GDP Growth and inflation from EUA 1948-2010
################################################################################

# import packages
library("readxl")
library("pracma") #for pseudo inverse pinv()
source("C:/Users/mathe/Desktop/Biblioteca/Trabalhos/R_Trabalhos/Bayesian Methods for Macroeconomists/FunctionsBMM.R")

# import data
data0 <- read_excel("C:/Users/mathe/Desktop/Biblioteca/Trabalhos/R_Trabalhos/Bayesian Methods for Macroeconomists/Code2017/CHAPTER2/DATA/DATAIN.XLS",
                    col_names = c("US_GDP_Growth", "Inflation"))
data = as.matrix(data0)

N = ncol(data)
horizon = 3
path = matrix(1, nrow = 3, ncol = 1) # constrained value for X in the forecast horizon
L = 2                                # VAR lag length
Y = data

# Take lags
X = c()
for (j in 1:L){
  X = cbind(X, lag0(data, j))
}
X = cbind(X, 1)
Y = Y[(L+1):nrow(Y),]
X = X[(L+1):nrow(X),]
T = nrow(X)

# OLS to be used in the IRF below
B = solve(t(X)%*%X)%*%t(X)%*%Y
res = Y - X%*%B
sigma = (t(res)%*%res)/T
A0 = chol(sigma)

# Calculate IRF to be used to construct R
# First shock
S = matrix(0, ncol = N)
S[1] = 1 # shock to the first equation
Z1 = irfsim(B,N,L,A0,S,horizon+L)
# Second shock
S = matrix(0, ncol = N)
S[2] = 1 # shock to the second equation equation
Z2 = irfsim(B,N,L,A0,S,horizon+L)

# Calculate unconditional forecast to be used to construct t
yhat1 = matrix(0, nrow = (horizon + L), ncol = N)
yhat1[(1:L),] = tail(Y, L)
for (i in (L+1):(horizon + L)){
  x = c()
  for (j in 1:L){
    x =  c(x, yhat1[(i-j),])
  }
  yhat1[i,] = c(x,1)%*%B
}
yhat1 = yhat1[(L+1):nrow(yhat1),]

# Construct the R matrix
R = rbind(c(Z1[1,2], Z2[1,2], 0, 0, 0, 0),
          c(Z1[2,2], Z2[2,2], Z1[1,2], Z2[1,2], 0, 0),
          c(Z1[3,2], Z2[3,2], Z1[2,2], Z2[2,2], Z1[1,2], Z2[1,2]))

# Construct the r matrix
r = path - yhat1[,2]

# Compute the restricted structural shocks
ehat = t(R)%*%pinv(R%*%t(R))%*%r
ehat = t(matrix(ehat, nrow = N, ncol = horizon))

# compute the conditional forecast
yhat2 = matrix(0, nrow = (horizon + L), ncol = N)
yhat2[(1:L),] = tail(Y, L)
for (i in (L+1):(horizon + L)){
  x = c()
  for (j in 1:L){
    x =  c(x, yhat2[(i-j),])
  }
  yhat2[i,] = c(x,1)%*%B + ehat[(i - L),]%*%A0
}
yhat2 = yhat2[(L+1):nrow(yhat2), ]

# Gibbs sampling algorithm
REPS = 5000; BURN = 3000
# Store forecast result 
forecastGDP = matrix(0, nrow = (REPS - BURN), ncol = (3))
forecastInflation = matrix(0, nrow = (REPS - BURN), ncol = (3))
yhatg = yhat2 # initialize conditional forecast
sig = sigma   # initialize error covariance

for(igibbs in 1:REPS){
  print(paste("Gibbs iteration: ", igibbs))
  
  # Step 1: Draw VAR parameters
  datag = rbind(data, yhatg) # append data
  YSTAR = datag
  
  # take lags
  XSTAR = c()
  for (j in 1:L){
    XSTAR = cbind(XSTAR, lag0(YSTAR, j))
  }
  XSTAR = cbind(XSTAR, 1)
  YSTAR = YSTAR[(L+1):nrow(YSTAR), ]
  XSTAR = XSTAR[(L+1):nrow(XSTAR), ]
  T = nrow(XSTAR)
  
  # Conditional mean
  M = matrix(solve(t(XSTAR)%*%XSTAR)%*%t(XSTAR)%*%YSTAR, ncol = 1)
  # Conditional variance
  V = kronecker(sig, solve(t(XSTAR)%*%XSTAR))
  
  bg = M + t(rnorm(N*((N*L) + 1))%*%chol(V))
  bg1 = matrix(bg, nrow = ((N*L )+ 1), ncol = N)
  
  # Draw sigma from the IW distribution
  e = YSTAR - XSTAR%*%bg1
  scale = t(e)%*%e
  sig = IWPQ(T, solve(scale))
  
  # A0 matrix
  A0g = chol(sig)
  
  #STep 2: Construct conditional forecast
  # impulse response
  # First shock
  S = matrix(0, ncol = N)
  S[1] = 1 # shock to the first equation
  Z1 = irfsim(bg1, N, L, A0, S, horizon+L)
  # Second shock
  S = matrix(0, ncol = N)
  S[2] = 1 # shock to the second equation equation
  Z2 = irfsim(bg1, N, L, A0, S, horizon+L)
  
  # Calculate unconditional forecast to be used to construct t
  yhat1 = matrix(0, nrow = (horizon + L), ncol = N)
  yhat1[(1:L),] = tail(Y, L)
  for (i in (L+1):(horizon + L)){
    x = c()
    for (j in 1:L){
      x =  c(x, yhat1[(i-j),])
    }
    yhat1[i,] = c(x,1)%*%bg1
  }
  yhat1 = yhat1[(L+1):nrow(yhat1),]
  
  # Construct the R matrix
  R = rbind(c(Z1[1,2], Z2[1,2], 0, 0, 0, 0),
            c(Z1[2,2], Z2[2,2], Z1[1,2], Z2[1,2], 0, 0),
            c(Z1[3,2], Z2[3,2], Z1[2,2], Z2[2,2], Z1[1,2], Z2[1,2]))
  
  # Construct the r matrix
  r = path - yhat1[,2]
  
  # Compute the mean of the distribution of restricted structural shocks
  MBAR = t(R)%*%pinv(R%*%t(R))%*%r
  
  # Compute the variance of the distribution of restricted structural shocks
  VBAR = t(R)%*%pinv(R%*%t(R))%*%R
  VBAR = diag(nrow = ncol(VBAR)) - VBAR
  
  # Draw structural shocks from the N(MBAR, VBAR) distribution
  edraw = MBAR + t(rnorm(nrow(MBAR))%*%Re(sqrtm(VBAR)$B))
  edraw = t(matrix(edraw, nrow = N, ncol = horizon))
  
  # Conditional forecast using new draw of shocks
  yhatg = matrix(0, nrow = (horizon + L), ncol = N)
  yhatg[1:L, ] = tail(Y, L)
  for (i in (L+1):(horizon + L)){
    x = c()
    for (j in 1:L){
      x =  c(x, yhatg[(i-j),])
    }
    yhatg[i,] = c(x,1)%*%bg1 + ehat[(i - L),]%*%A0g
  }
  yhatg = yhatg[(L+1):nrow(yhatg), ]
  
  if(igibbs > BURN){
    temp = igibbs - BURN
    forecastGDP[temp, ] = yhatg[ ,1]
    forecastInflation[temp, ] = yhatg[ ,2]
  }
}

# Plot of the forecast distribution
library(fanplot)
plot(data0$Inflation[1:210], type = "l", xlim = c(0, 270))
fan(data = forecastInflation, start = 210, p = seq(0.05, 0.95, 0.05), 
    fan.col = colorRampPalette(c("tomato", "gray90")))

plot(data0$US_GDP_Growth[1:210], type = "l", xlim = c(0, 270))
fan(data = forecastGDP, start = 210, p =seq(0.05, 0.95, 0.05), 
    fan.col = colorRampPalette(c("tomato", "gray90")))
