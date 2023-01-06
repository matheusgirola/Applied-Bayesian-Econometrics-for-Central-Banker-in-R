################################################################################
# VAR estimation via Gibbs sampling with steady state prior
# following Blake and Mumtaz, Chapter 2
# This is what they called "example3.m"
# Model to estimate: VAR(2) Y(t) = alpha + B(1)Y(t-1) + B(2)Y(t-2) + v(t)
# Y = quartely annual US GDP Growth and inflation from EUA 1948-2010

# import packages
library("readxl")
library("LaplacesDemon") #for Wishart distribution
source("C:/Users/MATHEUS/Desktop/Biblioteca/Trabalhos/R_Trabalhos/Bayesian Methods for Macroeconomists/FunctionsBMM.R")

# import data
data0 <- read_excel("C:/Users/MATHEUS/Desktop/Biblioteca/Trabalhos/R_Trabalhos/Bayesian Methods for Macroeconomists/Code2017/CHAPTER2/DATA/DATAIN.XLS",
                   col_names = c("US_GDP_Growth", "Inflation"))
data = as.matrix(data0)
Y = data

# Create Lag variables
N = ncol(Y)
L = 2 # lags
X = cbind(lag0(Y, 1), lag0(Y, 2), 1)
Y = Y[3:nrow(Y),]
X = X[3:nrow(X),]
T = nrow(X)

# Compute standard deviation of each series residual via an OLS Regression,
# to be used in setting the prior

# First variable
y = as.matrix(Y[,1], ncol = 1)
x = as.matrix(cbind(1 , X[,1]), ncol = 2)
b0 = solve(t(x)%*%x)%*%t(x)%*%y
s1 = sqrt((t(y - x%*%b0)%*%(y-x%*%b0))/(nrow(y)-2)) #std of residual standard error
s1 = s1[1] # return constante

# Second variable
y = as.matrix(Y[,2], ncol = 1)
x = as.matrix(cbind(1 , X[,2]), ncol = 2)
b0 = solve(t(x)%*%x)%*%t(x)%*%y
s2 = sqrt((t(y - x%*%b0)%*%(y-x%*%b0))/(nrow(y)-2)) #std of residual standard error
s2 = s2[1]

# Specify hyperparameters of the Minnesota priro
lamda1 = 1; lamda2 = 1; lamda3 = 1; lamda4 = 1

# Specify prior mean of the coefficients of the Two equations of the VAR
B01 = c(1,0,0,0)
B02 = c(0,1,0,0)
B0 = as.matrix(c(B01, B02), ncol = 1)

# Specify the prior variance of vec(B)
H = matrix(0, nrow = 8, ncol = 8)
# For equation 1 of the VAR
H[1,1] = (lamda1)^2                             # own lag
H[2,2] = ((s1*lamda1*lamda2)/s2)^2              # lag of other variable
H[3,3] = (lamda1/(2^lamda3))^2                  # own second lag
H[4,4] = ((s1*lamda1*lamda2)/(s2*(2^lamda3)))^2 # lag of other variable
# For equation 2 of the VAR 
H[5,5] = ((s2*lamda1*lamda2)/s1)^2              # lag of other variable
H[6,6] = (lamda1)^2                             # own lag
H[7,7] = ((s2*lamda1*lamda2)/(s1*(2^lamda3)))^2 # lag of other variable
H[8,8] = (lamda1/(2^lamda3))^2                  # own second lag

# Prior scale matrix for sigma the VAR covariance
S = diag(nrow = N) # S^bar

# Priro degrees of freedom
alpha = N + 1

# set priors for the long run mean which is a N by 1 vector
M0 = matrix(1, nrow = 1, ncol = 2)   # prior mean
V0 = diag(nrow = N)*0.001 # prior variance

#Starting values for gibbs sampling
# Starting values for mean with OLS
betaols = solve(t(X)%*%X)%*%t(X)%*%Y
F = rbind(t(betaols[1:(N*L),]), diag(1, nrow = N*(L-1), ncol = N*L)) # companion form
C = matrix(0, nrow = nrow(F), ncol = 1)                              
C[1:N] = t(betaols[((N*L) + 1),]) # vector of constant
MU = solve(diag(nrow = nrow(F)) - F)%*%C # Mean = (I - F)*C

e = Y - X%*%betaols
Sigma = (t(e)%*%e)/T

REPS = 10000; BURN = 5000

# to store forecast
forecastGDP = matrix(, nrow = (REPS - BURN), ncol = 42)
forecastInflation = matrix(, nrow = (REPS - BURN), ncol = 42)

for (i in 1:REPS){
  print(sprintf("%d gibbs step", i))
  # "demean" the data (take out mean)
  mean_matrix = c()
  for (k in 1:nrow(data)){
    mean_matrix = rbind(mean_matrix, t(MU[1:N]))
  }
  Y0 = data - mean_matrix #Y0 = Z - mean
  X0 = c()
  for (j in 1:L){
    X0 = cbind(X0, lag0(Y0, j))
  }
  Y0 = Y0[((L+1):nrow(Y0)),]
  X0 = X0[((L+1):nrow(X0)),]
  
  #step 1: draw the VAR coefficients
  bols = matrix(solve(t(X0)%*%X0)%*%t(X0)%*%Y0, ncol = 1)
  M = solve(solve(H) + kronecker(solve(Sigma), t(X0)%*%X0))%*%
    + (solve(H)%*%B0 + (kronecker(solve(Sigma), t(X0)%*%X0))%*%bols)
  V = solve(solve(H) + kronecker(solve(Sigma), t(X0)%*%X0))
  
  beta = M + t(rnorm(N*(N*L))%*%chol(V))
  beta1 = matrix(beta, nrow = (N*L), ncol = N)
  
  # Step 2: draw sigma from the Wishart distribution
  e = Y0 - X0%*%beta1
  # scale matrix
  scale = t(e)%*%e + S
  Sigma = IWPQ(T + alpha, solve(scale))#rinvwishart(nu = (T + alpha), S = solve(scale))
  
  # Step 3: draw MU, the long run mean 
  Y1 = Y - X[, 1:(ncol(X)-1)]%*%beta1 #Y= Zt - B1Zt-1 - ... - BpZt-p
  U = diag(nrow = N)
  jj = 1
  for (jx in 1:L){
    betai = beta1[jj:(jj + N -1),]
    U = rbind(U, t(betai))
    jj= jj + N
  }
  D = cbind(matrix(1, nrow = T, ncol = 1), -1*matrix(1, nrow = T, ncol = L))
  # posterior variance
  vstar1 = solve(t(U)%*%kronecker((t(D)%*%D), solve(Sigma))%*%U + solve(V0))
  # posterior mean
  mstar1 = vstar1%*%(t(U)%*%matrix((solve(Sigma)%*%t(Y1)%*%D), ncol = 1)+
                       + solve(V0)%*%t(M0))
  MU = mstar1 + t(rnorm(N)%*%chol(vstar1)) # draw Mu
  if (i > BURN){
    temp = i - BURN
    # forecast GDP Growth and Inflation for 3 years
    F = rbind(t(beta1[1:(N*L),]), diag(1, nrow = N*(L-1), ncol = N*L)) # companion form
    mu = c()
    for (i in 1:L){
      mu = rbind(mu, MU)
    }
    C = (diag(nrow = nrow(F)) - F)%*%mu
    yhat = matrix(0, nrow = 44, ncol = 2)
    yhat[1:2,] = tail(Y, n = 2)
    for (m in 3:44){
      yhat[m,] = t(C[1:N]) + c(yhat[(m-1),], yhat[(m-2),])%*%
        + matrix(beta, nrow = N*L, ncol = N) + rnorm(N)%*%chol(Sigma)
    }
    forecastGDP[temp,] = yhat[3:44, 1]
    forecastInflation[temp, ] = yhat[3:44, 2]
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
