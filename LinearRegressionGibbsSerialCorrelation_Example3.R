################################################################################
# Linear Regression estimation via Gibbs sampling, following Blake and Mumtaz 
# This is what they called "examplem3.m"
# Model to estimate: AR(2) Y(t) = alpha + B(1)Y(t-1) + B(2)y(t-2) + v(t)
#                          v(t) = rho*v(t-1) + epsilon(t)
# Y = quartely annual inflation from EUA 1948-2010

# import packages
library("readxl")
#library("tidyverse")
source("C:/Users/mathe/Desktop/Biblioteca/Trabalhos/R_Trabalhos/CodigoMumtaz/FunctionsMumtaz.R")

# import data
data <- read_excel("C:/Users/mathe/Desktop/Biblioteca/Trabalhos/R_Trabalhos/Bayesian Methods for Macroeconomists/Code2017/CHAPTER1/DATA/inflation.xls")
Y = as.matrix(data$Inflation)

#Create the lag variables
T = nrow(data)
Y_1 = lag0(Y, k = 1)
Y_2 = lag0(Y, k = 2) 

# Create Regressors
X = cbind(rep(1, T), Y_1, Y_2)

# remove missing obs
Y = as.matrix(Y[3:nrow(Y),])
X = as.matrix(X[3:nrow(X),])

# SET PRIORS AND STARTING VALUES

# Set priors for B
B_0 = c(0,0,0); Sigma_0 = diag(3)
# Set priors for sigma2
T_0 = 1; D_0 = 0.1
# Set priors for rho
rho0 = 0; Sigma0r = 1

# Set starting values]
B = B_0; rho = rho0; sigma2 = 1
REPS = 25000; BURNIN = 24000

# Variables to store outcomes
out_1 = matrix(, nrow = REPS, ncol = 3); out_2 = numeric(REPS); out_3 = numeric(REPS)
forecast = matrix(, nrow = 1000, ncol = 12)
B_sample = matrix(, nrow = 1000, ncol = 3)
sigma_sample =  numeric(1000)
rho_sample = numeric(1000)

# Gibbs sampling
for (i in 1:REPS){
  print(sprintf("%d gibbs step", i))
  # SAMPLE B CONDITIONAL ON RHO AND SIGMA
  # removal serial correlation
  ystar = Y - lag0(Y, k = 1)*rho
  xstar = X - lag0(X, k = 1)*rho
  ystar = ystar[2:nrow(ystar),]
  xstar = xstar[2:nrow(xstar),]
  
  M = solve(solve(Sigma_0) + (1/sigma2)*(t(xstar)%*%xstar))%*%(solve(Sigma_0)%*%B_0 + (1/sigma2)*(t(xstar)%*%ystar))
  V = solve(solve(Sigma_0) + (1/sigma2)*(t(xstar)%*%xstar))
  
  check = -1
  while (check < 0){            #check for stability
    B = M + t(rnorm(3)%*%chol(V))
    b = cbind(c(B[2], B[3]), c(1, 0))
    ee = max(abs(eigen(b)$values))
    if (ee <= 1){
      check = 1
    }
  }
  # COMPUTE RHO CONDITIONAL ON B AND SIGMA
  yv = Y - X%*%B
  xv = lag0(yv, k=1)
  yv = yv[2:length(yv),]
  xv = xv[2:length(xv),]
  
  MM = solve(solve(Sigma0r) + (1/sigma2)*(t(xv)%*%xv))%*%(solve(Sigma0r)%*%rho0 + (1/sigma2)*(t(xv)%*%yv))
  VV = solve(1/Sigma0r + (1/sigma2)*(t(xv)%*%xv))
  
  #draw rho and ensure stationarity
  check = -1
  while (check < 0){            #check for stability
    rho = MM + t(rnorm(1)%*%chol(VV))
    ee = abs(rho)
    if (ee <= 1){
      check = 1
    }
  }
  # store rho as scalar
  rho = rho[1]
  
  # SAMPLE SIGMA2 CONDITIONAL ON B AND RHO
  # Compute residuals
  resids = ystar - xstar%*%B
  # Compute posterior df and scale matrix
  T_1 = T_0 + T
  D_1 = D_0 + t(resids)%*%resids
  # Draw from IG
  z_0 = rnorm(T_1)
  z0z0 = t(z_0)%*%z_0
  sigma_2 = D_1/z0z0   #return a 1x1 array
  sigma_2 = sigma_2[1] #store it as a scalar
  
  out_1[i,] = t(B)
  out_2[i] = sigma_2
  out_3[i] = rho
  if (i > BURNIN){
    # add samples
    temp = i - BURNIN
    B_sample[temp,] = t(B)
    sigma_sample[temp] = sigma_2
    rho_sample[temp] = rho
    
    #Compute forecast for 2 years
    yhat = numeric(14)
    vhat = numeric(14)
    yhat[1:2] = tail(Y, n = 2)
    cfactor = sqrt(sigma_2)
    
    for (m in 3:14){
      vhat[m] = vhat[m-1]*rho + rnorm(1)*cfactor
      yhat[m] = c(1,yhat[m-1],yhat[m-2])%*%B+vhat[m]
    }
    forecast[temp,] = yhat[3:14]
  }
}

# Plot marginal posterior distribution
library("plotly")
fig1 <- plot_ly(x = B_sample[,1], type = "histogram")
fig2 <- plot_ly(x = B_sample[,2], type = "histogram")
fig3 <- plot_ly(x = B_sample[,3], type = "histogram")
fig4 <- plot_ly(x = sigma_sample, type = "histogram")
fig5 <- plot_ly(x = rho_sample, type = "histogram")

fig <- subplot(fig1, fig2, fig3, fig4, fig5, nrows = 3)
fig

# Check mean, standard deviations and credible intervals
mean = numeric(5)
confidence_int = matrix(, nrow = 5, ncol = 2)

j=1
for (i in list(B_sample[,1], B_sample[,2], B_sample[,3], sigma_sample, rho_sample)){
  model <- lm (i~1)
  mean[j] = model$coefficients
  confidence_int[j,] = confint(model, level = 0.90)
  j = j + 1
}
std = c(sd(B_sample[,1]),sd(B_sample[,2]),
        sd(B_sample[,3]),sd(sigma_sample), sd(rho_sample))

# Summary of parametes statistics
summary = as.data.frame(cbind(mean, std, confidence_int))
colnames(summary) = c("Mean", "Standard Deviation", 
                      "5% percentile", "95% percentile")
rownames(summary) = c("alpha", "Beta 1", "Beta 2", "sigma", "rho")
summary

# Plot of the forecast distribution
library(fanplot)
plot(data$Inflation[1:250], type = "l", xlim = c(0, 270))
fan(data =  forecast, start = 240)

#Check convergence of gibbs draws
plot(x = B_sample[,1], type = "l")
plot(x = B_sample[,2], type = "l")
plot(x = B_sample[,3], type = "l")
plot(x = sigma_sample, type = "l")
plot(x = rho_sample, type = "l")

acf(x = B_sample[,1])
acf(x = B_sample[,2])
acf(x = B_sample[,3])
acf(x = sigma_sample)
acf(x = rho_sample)
