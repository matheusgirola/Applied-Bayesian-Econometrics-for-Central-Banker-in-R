################################################################################
# VAR estimation with sign restrictions via Gibbs sampling, 
# following Blake and Mumtaz, Chapter 2
# The code is similar to example5, the difference is the algorithm to search the
# sign restriction matrix A0 is more efficient
# This is what they called "example6.m"
# Model to estimate: VAR(2) Y(t) = alpha + B(1)Y(t-1) + B(2)Y(t-2) + v(t)
# Y = (1) Federal Funds Rate; (2) Annual GDP growth; (3) Annual CPI Inflation;
# (4) Annual real consumption growth; (5) Unemployment rate; 
# (6) change in private investment; (7) net exports; (8) annual growth in M2;
# (9) 10 year government bond yield; (10) annual growth in stock prices;
# (11) annual growth in the yen dollar exchange rate, from EUA 1971-2010
################################################################################


# import packages
library("readxl")
library("LaplacesDemon") #for Wishart distribution

source("C:/Users/mathe/Desktop/Biblioteca/Trabalhos/R_Trabalhos/Bayesian Methods for Macroeconomists/FunctionsBMM.R")
REPS = 5000; BURN = 3000

# import data
data0 <- read_excel("C:/Users/mathe/Desktop/Biblioteca/Trabalhos/R_Trabalhos/Bayesian Methods for Macroeconomists/Code2017/CHAPTER2/DATA/usdata1.xlsx",
                    col_names = c("Fed Funds Rates", "GDP Growth", "Inflation",
                                  "Real Consumption Growth", "Unemployment",
                                  "Change in Private Investment", "Net Exports",
                                  "Growth in M2", "10 year Government bond",
                                  "Growth in stock prices", "Growth in yen-dollar"))
data = as.matrix(data0)
Y = data
N = ncol(Y)

#Create Lag variables
L = 2 # VAR lag length
X = c()
for (l in 1:L){
  X = cbind(X, lag0(Y, k = l))
}
X = cbind(X, 1)
Y = Y[(L+1):nrow(Y),]
X = X[(L+1):nrow(X),]
T = nrow(X)

# Prior for VAR coefficients
lamdaP = 1       # Tightness of prior on first lag
tauP = lamdaP*10 # Tightness of prior on the sum of coefficients
epsilonP = 1     # Tightness of the prior on the constant
muP = t(colMeans(Y))
sigmaP = numeric(N)
deltaP = numeric(N)
for (i in 1:N){
  ytemp = matrix(Y[,i], ncol = 1)
  xtemp = cbind(lag0(ytemp, 1), 1)
  ytemp = ytemp[2:nrow(ytemp),]
  xtemp = xtemp[2:nrow(xtemp),]
  btemp = solve(t(xtemp)%*%xtemp)%*%t(xtemp)%*%ytemp
  etemp = ytemp - xtemp%*%btemp
  stemp = (t(etemp)%*%etemp)/length(ytemp)
  deltaP[i] = btemp[1]
  sigmaP[i] = stemp
}

# Dummy data to implement priors
dummies = create_dummies(lamdaP, tauP, deltaP, epsilonP, L, muP, sigmaP, N)
yd = dummies$y
xd = dummies$x

# Append dummies to data
Y0 = rbind(Y, yd)
X0 = rbind(X, xd)

# Conditional Mean of the VAR coefficients
mstar = matrix(solve(t(X0)%*%X0)%*%t(X0)%*%Y0, ncol = 1) # OLS on the appended data
xx = t(X0)%*%X0
ixx = solve(t(xx)%*%xx)%*%t(xx)%*%diag(ncol(xx)) # inv(t(X0)X0) to be used later on Gibbs
sigma = diag(N) # starting values for sigma
out = array(0, dim = c((REPS-BURN),36, N)) # to store outcomes

for (i in 1:REPS){
  print(paste("Gibbs iteration: ", i))
  
  # Draw VAR coefficients
  vstar = kronecker(sigma, ixx)
  beta = mstar + t(rnorm(N*((N*L) + 1))%*%chol(vstar))
  
  # Draw covariance
  e = Y0 - X0%*%matrix(beta, nrow = ((N*L) + 1), ncol = N)
  scale = t(e)%*%e
  sigma = IWPQ((T +  nrow(yd)), solve(scale)) #rinvwishart(nu = (T +  nrow(yd)), S = solve(scale))
  
  if (i > BURN){
    temp = i - BURN
    
    # Impose sign restrictions
    chck = -1
    while (chck < 0){
      K = matrix(rnorm(N*N), ncol = N, nrow = N)
      # take QR decomposition of K. Note that Q is orthonormal
      Q = getqr(K) 
      # take the cholesnky decomposition of current sigma draw
      A0hat = chol(sigma)
      # Calculate the candidate matrix draw A0
      A0hat1 = (Q%*%A0hat)
      
      for(m in 1:N){
        #check signs
        e1 = A0hat1[m,1] > 0 # Response of R
        e2 = A0hat1[m,2] < 0 # Response of Y
        e3 = A0hat1[m,3] < 0 # Response of Inflation
        e4 = A0hat1[m,4] < 0 # Response of Consumption
        e5 = A0hat1[m,5] > 0 # Response of U
        e6 = A0hat1[m,6] < 0 # Response of investment
        e7 = A0hat1[m,8] < 0 # Response of money
        if ((e1 + e2 + e3 + e4 + e5 + e6 + e7)==7){
          MP = A0hat1[m,]
          chck = 10
        }else{
          # check signs but reverse them
          e1 = -A0hat1[m,1] > 0 # Response of R
          e2 = -A0hat1[m,2] < 0 # Response of Y
          e3 = -A0hat1[m,3] < 0 # Response of Inflation
          e4 = -A0hat1[m,4] < 0 # Response of Consumption
          e5 = -A0hat1[m,5] > 0 # Response of U
          e6 = -A0hat1[m,6] < 0 # Response of investment
          e7 = -A0hat1[m,8] < 0 # Response of money
          if ((e1 + e2 + e3 + e4 + e5 + e6 + e7)==7){
            MP = -A0hat1[m,]
            chck = 10
          }
        }
      }
    }
    # re-shuffle rows of A0hat1 and insert MP in the first row
    A0x = c() # will hold rows of A0hat1 not equal to MP
    for (m in 1:N){
      ee = sum(abs(A0hat1[m,]) == abs(MP))
      if (ee == 0){
        A0x = rbind(A0x, A0hat1[m,])
      }
    }
    A0new = rbind(MP, A0x) # A0 to be used in IRF
    yhat = matrix(0, nrow = 36, ncol = N)
    vhat = matrix(0, nrow = 36, ncol = N)
    vhat[3,1] = 1 # Shock to the Fed fund rates
    for (j in 3:36){
      yhat[j,] = c(yhat[(j-1),], yhat[(j-2),], 0)%*%
        + matrix(beta, nrow=((N*L)+1), ncol = N) +vhat[j,]%*%A0new
    }
    out[temp, , ] = yhat
  }
}

# PLOT IRFS ####################################################################
plotIRF <- function(col_data){
  summary = matrix(0, nrow = 3, ncol = 36)
  rownames(summary) = c("16% percentile", "Median", "68% percentile")
  for (m in 1:36){
    summary[,m] = quantile(out[,m,col_data], c(0.16, 0.50, 0.68))
  }
  summary = rbind(summary, 0)
  matplot(x = c(1, 36), y = c(min(summary), max(summary)), type =  "n", 
          xlab = "Months after schock", ylab = "Response",
          main = paste("IRF", names(data0)[col_data] ))
  matpoints(x = (1:36), y = t(summary), type = "l", lty = 1, lwd = 1.5,)
  legend(x = "bottomright", c("16% percentile", "Median", "68% percentile", "0 line"),
         fill = c("green", "red", "black", "blue"), cex = 0.5)
}

plotIRF(1)  # Fed Fund Rates
plotIRF(2)  # GDP Growht
plotIRF(3)  # Inflation
plotIRF(4)  # Real comsumption Growth
plotIRF(5)  # Unemployment
plotIRF(6)  # Change in Private Investment
plotIRF(7)  # Net Export
plotIRF(8)  # Growth in M2
plotIRF(9)  # 10 year Governemnt Bond
plotIRF(10) # Growth in stock prices
plotIRF(11) # Growth in yen-dollar exchange