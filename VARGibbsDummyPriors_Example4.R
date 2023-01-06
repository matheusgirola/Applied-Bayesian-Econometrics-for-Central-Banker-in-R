################################################################################
# VAR estimation with dummy priors via Gibbs sampling, 
# following Blake and Mumtaz, Chapter 2
# This is what they called "example4.m"
# Model to estimate: VAR(2) Y(t) = alpha + B(1)Y(t-1) + B(2)Y(t-2) + v(t)
# Y = quartely annual US GDP Growth and inflation from EUA 1948-2010

# import packages
library("readxl")
library("LaplacesDemon") #for Wishart distribution
source("C:/Users/mathe/Desktop/Biblioteca/Trabalhos/R_Trabalhos/CodigoMumtaz/FunctionsMumtaz.R")

# import data
data <- read_excel("C:/Users/mathe/Desktop/Biblioteca/Trabalhos/R_Trabalhos/Bayesian Methods for Macroeconomists/Code2017/CHAPTER2/DATA/DATAIN.XLS",
                   col_names = c("US_GDP_Growth", "Inflation"))
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

# mean of the data
mu = colMeans(Y)

# specify parameters of the minessota prior
tau = 0.1  # controls prior on own first lag
d = 1      # decay for higher lags
lamdac = 1 # tightness for prior on constant
lamda = 1  # tightness on prior of sum of coefficients
delta = 1  # tightness of cointegration prior

# Specify dummies observations for first lag
yd1 = matrix(c((1/tau)*s1, 0, 0, (1/tau)*s2), nrow = 2, ncol = 2)
xd1 = rbind(c(0, (1/tau)*s1, 0, 0, 0), c(0, 0, (1/tau)*s2, 0, 0))

# Specify dummies for second lag
yd2 = matrix(0, ncol = 2, nrow = 2)
xd2 = rbind(c(0, 0, 0, (1/tau)*s1*2^d, 0), c(0, 0, 0, 0, (1/tau)*s2*2^d))

# Specify dummies for the constants
yd3 = matrix(0, ncol = 2, nrow = 2)
xd3 = rbind(c(1/lamdac, 0, 0, 0, 0), c(1/lamdac, 0, 0, 0, 0))

# Specify dummies for the sum of coefficients
yd4 = rbind(c(lamda*mu[1], 0), c(0, lamda*mu[2]))
xd4 = rbind(c(0, lamda*mu[1], 0, lamda*mu[1], 0),
            c(0, 0, lamda*mu[2], 0, lamda*mu[2]))

# Specify common stochastic trend dummies
yd5 = matrix(c(delta*mu[1], delta*mu[2]), nrow = 1)
xd5 = matrix(c(delta , delta*mu[1], delta*mu[2], delta*mu[1], delta*mu[2]), nrow =1)

# Specify dummies for covariance matrix
yd6 = rbind(c(s1, 0), c(0, s2))
xd6 = matrix(0, nrow = 2, ncol = 5)

# All dummy observations
yd = rbind(yd1, yd2, yd3, yd4, yd5, yd6)
xd = rbind(xd1, xd2, xd3, xd4, xd5, xd6)

# Append dummies to data
Ystar = rbind(Y, yd)
Xstar = rbind(X, xd)
Tstar = nrow(Xstar)

# Compute posterior mean
betahat = solve(t(Xstar)%*%Xstar)%*%t(Xstar)%*%Ystar

# Compute initial value of sigma
e = Ystar - Xstar%*%betahat
sigma = (t(e)%*%e)/Tstar

# Matrices to store forecast
REPS = 2000; BURN = 1000
forecastGDP = matrix(0, nrow = (REPS - BURN), ncol = 12)
forecastInflation = matrix(0, nrow = (REPS - BURN), ncol = 12)

for (i in 1:REPS){
  print(sprintf("%d gibbs step", i))
  
  M = matrix(betahat, ncol = 1)
  V = kronecker(sigma, solve(t(Xstar)%*%Xstar))
  
  # Draw beta
  beta = M + t(rnorm(N*((N*L) + 1))%*%chol(V))
  
  # Draw Sigma
  e = Ystar - Xstar%*%matrix(beta, nrow = ((N*L) + 1), ncol = N)
  scale = t(e)%*%e
  sigma = rinvwishart(nu = (Tstar), S = solve(scale))
  
  if (i > BURN){
    temp = i - BURN
    yhat = matrix(0, nrow = 14, ncol = 2)
    yhat[1:2,] = tail(Y, n =2)
    for (j in 3:14){
      yhat[j,] = c(1, yhat[j - 1,], yhat[j - 2,])%*%
        matrix(beta, nrow = ((N*L) + 1), ncol = N) + rnorm(N)%*%chol(sigma)
    }
    forecastGDP[temp,] = yhat[3:14, 1]
    forecastInflation[temp,] = yhat[3:14, 2]
  }
}

# Plot of the forecast distribution
library(fanplot)
plot(data$Inflation[0:250], type = "l", xlim = c(0, 270))
fan(data = forecastInflation, start = 250, p =seq(0.1, 0.9, 0.1))

library(fanplot)
plot(data$US_GDP_Growth[0:250], type = "l", xlim = c(0, 270))
fan(data = forecastGDP, start = 250, p =seq(0.1, 0.9, 0.1))