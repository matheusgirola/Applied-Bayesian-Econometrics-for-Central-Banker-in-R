################################################################################
# VAR estimation via Gibbs sampling, following Blake and Mumtaz, Chapter 2
# This is what they called "example2.m"
# Model to estimate: VAR(2) Y(t) = alpha + B(1)Y(t-1) + B(2)Y(t-2) + v(t)
# Y = federal fund rates, 10 year government bond yield, unemployment rate
# and annual CPI inflation from EUA 2007 - 2010

# import packages
library("readxl")
library("LaplacesDemon") #for Wishart distribution
source("C:/Users/mathe/Desktop/Biblioteca/Trabalhos/R_Trabalhos/CodigoMumtaz/FunctionsMumtaz.R")

# import data
data <- read_excel("C:/Users/mathe/Desktop/Biblioteca/Trabalhos/R_Trabalhos/Bayesian Methods for Macroeconomists/Code2017/CHAPTER2/DATA/dataUS.XLS",
                   col_names = c("Federal Fund Rates", "10 Year Government Bond Yield",
                                "Unemployment Rate","Inflation"))
Y = as.matrix(data)

# Create Lag variables
N = ncol(Y)
L = 2 # lags
X = cbind(1, lag0(Y, 1), lag0(Y, 2))
Y = Y[3:nrow(Y),]
X = X[3:nrow(X),]
T = nrow(X)

# Compute standard deviation of each series residual via an OLS Regression,
# to be used in setting the prior

# First variable
y = as.matrix(Y[,1], ncol = 1)
x = X[,1:2]
b0 = solve(t(x)%*%x)%*%t(x)%*%y
s1 = sqrt((t(y - x%*%b0)%*%(y-x%*%b0))/(nrow(y)-2)) #std of residual standard error
s1 = s1[1] # return constante

# Second variable
y = as.matrix(Y[,2], ncol = 1)
x = X[,c(1,3)]
b0 = solve(t(x)%*%x)%*%t(x)%*%y
s2 = sqrt((t(y - x%*%b0)%*%(y-x%*%b0))/(nrow(y)-2)) #std of residual standard error
s2 = s2[1]

# Third variable
y = as.matrix(Y[,3], ncol = 1)
x = X[,c(1,4)]
b0 = solve(t(x)%*%x)%*%t(x)%*%y
s3 = sqrt((t(y - x%*%b0)%*%(y-x%*%b0))/(nrow(y)-2)) #std of residual standard error
s3 = s3[1]

# Fourth variable
y = as.matrix(Y[,4], ncol = 1)
x = X[,c(1,5)]
b0 = solve(t(x)%*%x)%*%t(x)%*%y
s4 = sqrt((t(y - x%*%b0)%*%(y-x%*%b0))/(nrow(y)-2)) #std of residual standard error
s4 = s4[1]

# Parameters to control the prior
lamda1 = 0.1  #Tightness of prior on AR coefficients
lamda3 = 0.05 #Tightness of prior in higher lags
lamda4 = 1    #Tightness prior on the constante term

# Specify the prior mean of the coefficients on the two VAR equations
B0 = matrix(0, nrow = ((N*L)+1), ncol = (N))
for (i in 1:N){
  B0[i + 1, i] = 0.95
}
B0 = matrix(B0, ncol = 1) #vec(B)

# Specify the prior variance of vec(B)
H = diag(nrow = (N*((N*L) + 1)))
# small for coefficients we want close to zero
H[3,3] = 1e-09 #for b12
H[4,4] = 1e-09 #for b13
H[5,5] = 1e-09 #for b14
H[7,7] = 1e-09 #for d12
H[8,8] = 1e-09 #for d13
H[9,9] = 1e-09 #for d14
# for the other coefficients we use normal conjugate prior
# first equation
H[1,1] = (s1*lamda4)^2
H[2,2] = (lamda1)^2
H[6,6] = (lamda1/(2^lamda3))^2
# second equation
H[10,10] = (s2*lamda4)^2
H[11,11] = ((s2*lamda1)/s1)^2
H[12,12] = (lamda1)^2
H[13,13] = ((s2*lamda1)/s3)^2
H[14,14] = ((s2*lamda1)/s4)^2
H[15,15] = ((s2*lamda1)/(s1*(2^lamda3)))^2
H[16,16] = (lamda1/(2^lamda3))^2
H[17,17] = ((s2*lamda1)/(s3*(2^lamda3)))^2
H[18,18] = ((s2*lamda1)/(s4*(2^lamda3)))^2
# third equation
H[19,19] = (s3*lamda4)^2
H[20,20] = ((s3*lamda1)/s1)^2
H[21,21] = ((s3*lamda1)/s2)^2
H[22,22] = (lamda1)^2
H[23,23] = ((s3*lamda1)/s4)^2
H[24,24] = ((s3*lamda1)/(s1*(2^lamda3)))^2
H[25,25] = ((s3*lamda1)/(s2*(2^lamda3)))^2
H[26,26] = (lamda1/(2^lamda3))^2
H[27,27] = ((s3*lamda1)/(s4*(2^lamda3)))^2
# fourth equation
H[28,28] = (s4*lamda4)^2
H[29,29] = ((s4*lamda1)/s1)^2
H[30,30] = ((s4*lamda1)/s2)^2
H[31,31] = ((s4*lamda1)/s3)^2
H[32,32] = (lamda1)^2
H[33,33] = ((s4*lamda1)/(s1*(2^lamda3)))^2
H[34,34] = ((s4*lamda1)/(s2*(2^lamda3)))^2
H[35,35] = ((s4*lamda1)/(s3*(2^lamda3)))^2
H[36,36] = (lamda1/(2^lamda3))^2

# Prior scale matrix for sigma the VAR covariance
S = diag(nrow = N)

# Prior degrees of freedom
alpha = N + 1

# Starting values for the gibbs sampling
Sigma = diag(nrow = N)
betaols = (solve(t(X)%*%X)%*%t(X)%*%Y)
betaols = matrix(betaols, nrow = length(betaols), ncol = 1)
REPS = 40000
BURN = 30000
IRF_R  = matrix(0, nrow = (REPS - BURN), ncol = 58)
IRF_GB = matrix(0, nrow = (REPS - BURN), ncol = 58)
IRF_U  = matrix(0, nrow = (REPS - BURN), ncol = 58)
IRF_P  = matrix(0, nrow = (REPS - BURN), ncol = 58)

for (i in 1:REPS){
  print(sprintf("%d gibbs step", i))
  # Step 1: Draw the VAR coefficients
  M = solve(solve(H) + kronecker(solve(Sigma), t(X)%*%X))%*%
    + (solve(H)%*%B0 + (kronecker(solve(Sigma), t(X)%*%X))%*%betaols)
  V = solve(solve(H) + kronecker(solve(Sigma), t(X)%*%X))
  # check stability of the VAR
  check = -1
  while(check < 0){
    beta = M + t(rnorm(N*((N*L) + 1))%*%chol(V))
    CH = stabilityABE(beta,N,L)
    if (CH == 0){
      check = 1
    }
  }
  
  # Step 2: draw sigma from the IW distribution
  e = Y - X%*%matrix(beta, nrow = ((N*L) + 1), ncol = N)
  # scale matrix
  scale = t(e)%*%e + S
  Sigma = rinvwishart(nu = (T + alpha), S = solve(scale))
  if (i > BURN){
    temp = i - BURN
    # IRF using Cholensky Decomposition
    A0 = chol(Sigma)
    v = matrix(0, nrow = 60, ncol = N)
    v[(L+1), 2] = -1 #schock on government bondyield
    yhat = matrix(0, nrow = 60, ncol = N)
    for (j in 3:60){
      yhat[j,] = c(0, yhat[(j-1),], yhat[(j-2),])%*%
        + matrix(beta, nrow = ((N*L) + 1), ncol = N) + v[j,]%*%A0
    }
    IRF_R[temp,]  = yhat[3:60,1]
    IRF_GB[temp,]  = yhat[3:60,2]
    IRF_U[temp,]   = yhat[3:60,3]
    IRF_P[temp,]   = yhat[3:60,4]
  }
}

# Plot IRF of each variable
# Federal Fund Rates
summary = matrix(0, nrow = 3, ncol = 58)
rownames(summary) = c("16% percentile", "Median", "68% percentile")
for (m in 1:58){
  summary[,m] = quantile(IRF_R[,m], c(0.16, 0.50, 0.68))
}
summary = rbind(summary, 0)
matplot(x = c(1, 58), y = c(min(summary), max(summary)), type =  "n", 
        xlab = "Months after schock", ylab = "Response",
        main = paste("IRF", names(data)[1] ))
matpoints(x = (1:58), y = t(summary), type = "l", lty = 1, lwd = 1.5,)
legend(x = "bottomright", c("16% percentile", "Median", "68% percentile", "0 line"),
       fill = c("green", "red", "black", "blue"), cex = 0.5)

# 10 year government yield
summary = matrix(0, nrow = 3, ncol = 58)
rownames(summary) = c("16% percentile", "Median", "68% percentile")
for (m in 1:58){
  summary[,m] = quantile(IRF_GB[,m], c(0.16, 0.50, 0.68))
}
summary = rbind(summary, 0)
matplot(x = c(1, 58), y = c(min(summary), max(summary)), type =  "n", 
        xlab = "Months after schock", ylab = "Response",
        main = paste("IRF", names(data)[2] ))
matpoints(x = (1:58), y = t(summary), type = "l", lty = 1, lwd = 1.5,)
legend(x = "bottomright", c("16% percentile", "Median", "68% percentile", "0 line"),
       fill = c("green", "red", "black", "blue"), cex = 0.5)

# Unemployment
summary = matrix(0, nrow = 3, ncol = 58)
rownames(summary) = c("16% percentile", "Median", "68% percentile")
for (m in 1:58){
  summary[,m] = quantile(IRF_U[,m], c(0.16, 0.50, 0.68))
}
summary = rbind(summary, 0)
matplot(x = c(1, 58), y = c(min(summary), max(summary)), type =  "n", 
        xlab = "Months after schock", ylab = "Response",
        main = paste("IRF", names(data)[3] ))
matpoints(x = (1:58), y = t(summary), type = "l", lty = 1, lwd = 1.5,)
legend(x = "bottomright", c("16% percentile", "Median", "68% percentile", "0 line"),
       fill = c("green", "red", "black", "blue"), cex = 0.5)

# Inflation
summary = matrix(0, nrow = 3, ncol = 58)
rownames(summary) = c("16% percentile", "Median", "68% percentile")
for (m in 1:58){
  summary[,m] = quantile(IRF_P[,m], c(0.16, 0.50, 0.68))
}
summary = rbind(summary, 0)
matplot(x = c(1, 58), y = c(min(summary), max(summary)), type =  "n", 
        xlab = "Months after schock", ylab = "Response",
        main = paste("IRF", names(data)[4] ))
matpoints(x = (1:58), y = t(summary), type = "l", lty = 1, lwd = 1.5,)
legend(x = "bottomright", c("16% percentile", "Median", "68% percentile", "0 line"),
       fill = c("green", "red", "black", "blue"), cex = 0.5)