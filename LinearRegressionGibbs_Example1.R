################################################################################
# Linear Regression estimation via Gibbs sampling, following Blake and Mumtaz 
# This is what they called "examplem1.m"
# Model to estimate: AR(2) Y(t) = alpha + B(1)Y(t-1) + B(2)y(t-2) + v(t)
# Y = quartely annual inflation from EUA 1948-2010

# import packages
library("readxl")
#library("tidyverse")


# import data
data <- read_excel("C:/Users/mathe/Desktop/Biblioteca/Trabalhos/R_Trabalhos/Bayesian Methods for Macroeconomists/Code2017/CHAPTER1/DATA/inflation.xls")
Y = data$Inflation

#Create the lag variables
T = nrow(data)
Y_1 = append(diff(Y, lag = 1), 0, after = 0) 
Y_2 = append(diff(Y, lag = 2), c(0,0), after = 0)  

# Create Regressors
X = cbind(rep(1, T), Y_1, Y_2)

# SET PRIORS AND STARTING VALUES

# Set priors for B
B_0 = c(0,0,0); Sigma_0 = diag(3)

# Set priors for sigma2
T_0 = 1; D_0 = 0.1

# Starting Values
B = B_0; sigma_2 = 1
REPS = 10000; BURNIN = 9000
out_1 = matrix(, nrow = REPS, ncol = 3); out_2 = numeric(REPS)
t_X = t(X)

for (i in 1:REPS){
# SAMPLE B CONDITIONAL ON SIGMA(M*,V*)
M = solve(solve(Sigma_0) + (1/sigma_2)*(t_X%*%X))%*%(solve(Sigma_0)%*%B_0 + (1/sigma_2)*(t_X%*%Y))
V = solve(solve(Sigma_0) + (1/sigma_2)*t_X%*%X)

check = -1
while (check < 0){            #check for stability
                  B = M + t(rnorm(3)%*%chol(V))
                  b = cbind(c(B[2], B[3]), c(1, 0))
                  ee = max(abs(eigen(b)$values))
                  if (ee <= 1){
                              check = 1
                  }
}

# Sample sigma 2 conditional on B from IG(T_1,D_1)
# Compute residuals
resids = Y - X%*%B
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

sprintf("%d gibbs step", i)
}

# Remove Burnin
B_samples = out_1[-1:-BURNIN,]
sigma_samples = out_2[-1:-BURNIN]

# Plot marginal posterior distribution
library("plotly")
fig1 <- plot_ly(x = B_samples[,1], type = "histogram")
fig2 <- plot_ly(x = B_samples[,2], type = "histogram")
fig3 <- plot_ly(x = B_samples[,3], type = "histogram")
fig4 <- plot_ly(x = sigma_samples, type = "histogram")

fig <- subplot(fig1, fig2, fig3, fig4, nrows = 2)
fig


# Check mean, standard deviations and credible intervals

mean = numeric(4)
confidence_int = matrix(, nrow = 4, ncol = 2)

j=1
for (i in list(B_samples[,1], B_samples[,2], B_samples[,3], sigma_samples)){
  model <- lm (i~1)
  mean[j] = model$coefficients
  confidence_int[j,] = confint(model, level = 0.90)
  j = j + 1
}

std = c(sd(B_samples[,1]),sd(B_samples[,2]),
        sd(B_samples[,3]),sd(sigma_samples))

summary = as.data.frame(cbind(mean, std, confidence_int))
colnames(summary) = c("Mean", "Standard Deviation", 
                      "5% percentile", "95% percentile")
summary
