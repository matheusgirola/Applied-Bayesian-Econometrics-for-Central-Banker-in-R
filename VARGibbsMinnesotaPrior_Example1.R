################################################################################
# VAR estimation via Gibbs sampling, following Blake and Mumtaz, Chapter 2
# This is what they called "example1.m"
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

# Specify hyperparameters of the Minnesota priro
lamda1 = 1; lamda2 = 1; lamda3 = 1; lamda4 = 1

# Specify prior mean of the coefficients of the Two equations of the VAR
B01 = c(0,1,0,0,0)
B02 = c(0,0,1,0,0)
B0 = as.matrix(c(B01, B02), ncol = 1)

# Specify the prior variance of vec(B)
H = matrix(0, nrow = 10, ncol =10)
# For equation 1 of the VAR
H[1,1] = (s1*lamda4)^2                          # constant
H[2,2] = (lamda1)^2                             # own lag
H[3,3] = ((s1*lamda1*lamda2)/s2)^2              # lag of other variable
H[4,4] = (lamda1/(2^lamda3))^2                  # own second lag
H[5,5] = ((s1*lamda1*lamda2)/(s2*(2^lamda3)))^2 # lag of other variable
# For equation 2 of the VAR 
H[6,6] = (s2*lamda4)^2                          # constant
H[7,7] = ((s2*lamda1*lamda2)/s1)^2              # lag of other variable
H[8,8] = (lamda1)^2                             # own lag
H[9,9] = ((s2*lamda1*lamda2)/(s1*(2^lamda3)))^2 # lag of other variable
H[10,10] = (lamda1/(2^lamda3))^2                # own second lag

# Prior scale matrix for sigma the VAR covariance
S = diag(nrow = N) # S^bar

# Priro degrees of freedom
alpha = N + 1

#Starting values for gibbs sampling
Sigma = diag(nrow = N)
betaols = matrix((solve(t(X)%*%X)%*%t(X)%*%Y), nrow = 10, ncol = 1)
REPS = 10000
BURN = 9000
#forecastGDP = matrix(0, nrow = 1000, ncol = 12)
#forecastInflation = matrix(0, nrow = 1000, ncol = 12)
forecastGDP = matrix(0, ncol = 1000, nrow = (nrow(Y) + 12))
forecastInflation = matrix(0, ncol = 1000, nrow = (nrow(Y) + 12))
for (i in 1:REPS){
  print(sprintf("%d gibbs step", i))
  # Step 1: Draw the VAR coefficients
  M = solve(solve(H) + kronecker(solve(Sigma), t(X)%*%X))%*%
    + (solve(H)%*%B0 + (kronecker(solve(Sigma), t(X)%*%X))%*%betaols)
  V = solve(solve(H) + kronecker(solve(Sigma), t(X)%*%X))
  beta = M + t(rnorm(N*((N*L) + 1))%*%chol(V))
  
  # Draw Sigma from the I distribution
  e = Y - X%*%matrix(beta, nrow = ((N*L) + 1), ncol = N)
  # scale matrix
  scale = (t(e)%*%e) + S
  Sigma = rinvwishart(nu = (T + alpha), S = solve(scale))
  
  if ( i > BURN){
    temp = i - BURN
    # Forecast GDP Growth and Inflation
    yhat = matrix(0, nrow = 14, ncol =2)
    yhat[1:2,] = tail(Y, n = 2)
    for (m in 3:14){
      
      yhat[m,] = c(1, yhat[(m - 1),], yhat[(m -2),])%*%
        + matrix(beta, nrow = ((N*L)+1), ncol = N) + rnorm(N)%*%chol(Sigma)
    }
    forecastGDP[,temp] = c(Y[,1], yhat[3:14, 1])
    forecastInflation[,temp] = c(Y[,2], yhat[3:14, 2])
  }
}

# Plot of the forecast distribution
library(fanplot)
plot(data$Inflation[200:240], type = "l", xlim = c(0, 60))
fan(data = forecastInflation, start = 240, p =seq(0.1, 0.9, 0.1), type = "interval")

plot(data$US_GDP_Growth[1:240], type = "l", xlim = c(0, 270))
fan(data = forecastGDP, start = 240, p =seq(0.05, 0.95, 0.05), 
    fan.col = colorRampPalette(c("tomato", "gray90")))