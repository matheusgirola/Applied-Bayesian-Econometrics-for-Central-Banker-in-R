################################################################################
# Kalman Filter for State Space model with time varying coefficients
# following Blake and Mumtaz, Chapter 3
# 
# This is what they called "example1.m"
# 
# Model to estimate: Y(t) = Beta(t)*X(t) + e1
#                    Beta(t) = mu + F*Beta(t-1) + e2
# Where the data is artificially generated
################################################################################

# Generate artifical data
t = 500 
Q = 0.001
R = 0.01
F = 1; mu =0 # These are fixed
e1 = rnorm(t)*sqrt(R)
e2 = rnorm(t)*sqrt(Q)
Beta = matrix(0, nrow = t)
Y = matrix(0, nrow = t)
X = rnorm(t)
for (j in 2:t){
  Beta[j, ] = Beta[(j-1), ] + e2[j]       # transition equation
  Y[j, ] =  X[j]%*%t(Beta[j, ]) + e1[j] # observation equation
}

# Start ofthe Kalman Filter ####################################################
# Step 1: set up matrices for the Kalman Filter
beta00 = matrix(0)               # State variable b[0|0]
p00 = 1                          # variance of the state variable p[0|0]
beta_tt = c()                    # will hold the filtered state variable
ptt = array(0, dim = c(t, 1, 1)) # will hold its variance

# Initialize state variable
beta11 = beta00 # b[t-1|t-1]
p11 = p00       # p[t-1|t-1]

#Loop from period 1 to end of the sample
for (i in 1:t){
  x = X[i]
  
  # Look at equations 3.9 on the book
  # Prediction
  beta10 = mu + beta11%*%F
  p10 = F%*%p11%*%t(F) + Q
  yhat = t(x%*%t(beta10))
  eta = Y[i, ] - yhat
  feta = (x%*%p10%*%t(x)) + R
  
  # Updating
  K = (p10%*%t(x))%*%solve(feta) # Kalman gain
  beta11 = t(t(beta10) + K%*%t(eta))
  p11 = p10 - K%*%(x%*%p10)
  
  ptt[i, , ] =  p11
  beta_tt = rbind(beta_tt, beta11)
}

# Plot both the artificial data and the estimated one using Kalman Filter
g <- ggplot()
g <- g + geom_line(aes(x = (1:500), y = Beta, colour = "Data"))
g <- g + geom_line(aes(x = (1:500), y = beta_tt, colour = "Estimated"))
g <- g + xlab("Period t") + ylab("Coefficients") + theme(legend.title = element_blank())
g


