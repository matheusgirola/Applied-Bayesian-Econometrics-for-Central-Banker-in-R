#carterkohnvar##################################################################
carterkohnvar <- function(Y,X,Q,iamat,hlast,beta0,P00,L,CHECK,maxdraws,EX){
  T = nrow(Y)
  N = ncol(Y)
  
  #step 2a: Set up matrices for the kalman filter
  ns = ncol(beta0)
  mu = 0
  beta11 = beta0
  p11 = P00
  
  #step 2b: run kalman filter
  for (i in 1:T){
    x = t(kronecker(diag(1,N), X[i,]))
    H = diag(hlast[i+1,])
    R = iamat%*%H%*%t(iamat)
    
    # Prediction
    beta10 = beta11
    p10 = p11
    yhat = t(x%*%t(beta10))
    eta = Y[i,] - yhat
    feta =(x%*%p10%*%t(x)) + R
    # Updating
    K = (p10%*%t(x))%*%invpd(feta)
    beta11 = t(t(beta10) + K%*%t(eta))
    p11 = p10 - K%*%(x%*%p10)
    ptt = p11
    beta_tt = beta11
  }
  rm(list = c("beta11", "p11"))
  # End of Kalman Filter

  # step 2c: Backward recursion to calculate the mean and variance of the
  # state vector distribution
  chck = -1
  problem = 0
  trys = 1
  
  while ((chck<0) && (trys <=maxdraws)){
    wa = rnorm(ns)
    beta2 = beta_tt + (wa%*%cholx(ptt)) #draw for beta in period t from N(beta_tt, ptt)
    error = Y - X%*%matrix(beta2, nrow=((N*L)+EX), ncol=N) #var residuals
    
    roots = stability(t(beta2), N, L, EX)
    
    if (CHECK){
      if (sum(roots)==0){
        chck=1
      }
      else{
        trys=trys+1
      }
    }
    else{
      chck =1
    }
  }
  if (CHECK){
    if (chck<0){
      problem=1
    }
  }
  
  return(list(beta2 = beta2, error = error, roots = roots, problem = problem))
}
#chofac#########################################################################
chofac <- function(N,chovec){
  # Written by Tim Cogley
  # This transforms a vector of cholensky coefficients - chovec - into a lower
  # triangular cholesnky matrix - CF - of size N
  
  CF = diag(1, N)
  i = 1
  for (j in 2:N){
    k = 1
    while (k < j){
      CF[j,k] = chovec[i,1]
      i = i + 1
      k = k +1
    }
  }
  return(CF)
}

#cholx##########################################################################
cholx <- function(x){
  #library("expm") # for sqrtm() function
  # check if cholesnky factorization occurs succesfully. If not, returns the
  # transpose of the real parte of the matrix square root of x
  step1 =sqrtm(x)$B
  #print(typeof(step1))
  step2 =Re(step1)
  out = t(step2)
  #out = t(Re(sqrtm(x)))
  
  return(out)
}
#create_dummies#################################################################
create_dummies <- function (lamda,tau,delta,epsilon,p,mu,sig,n,ph,epsilonH){
  # Creates matrices of dummy observations, following Banbura et al. 2007
  # lamda tightness parameter
  # tau  prior on sum of coefficients
  # delta prior mean for VAR coefficients
  # epsilon tigtness of the prior around constant
  # mu sample mean of the data
  # sig vector of AR residual variances for the data, input from getdummies()
  
  #Initialize vectors
  yd1 = c()
  yd2 = c()
  xd1 = c()
  xd2 = c()
  
  #Get dummy matrices in equation (5) of Banbura et al. 2007
  #yd
  if (lamda > 0){
    if (epsilon > 0){
      yd1 = rbind(diag(sig*delta)/lamda, 
                  matrix(0,nrow =(n*(p-1)), ncol = n))
      yd1 = rbind(yd1, diag(sig))
      yd1 = rbind(yd1, 0)
      
      #xd
      jp = diag(1:p)
      xd1 = cbind(kronecker(jp,(diag(sig)/lamda)), 0)
      xd1 = rbind(xd1, 
                  matrix(0, nrow = n, ncol = ((n*p)+1)),
                  cbind(matrix(0, nrow = 1, ncol = (n*p)), epsilon))
    }else{
      #yd
      yd1 = rbind(diag(sig*delta)/lamda, 
                  matrix(0,nrow =(n*(p-1)), ncol = n))
      yd1 = rbind(yd1, diag(sig))
      
      #xd
      jp = diag(1:p)
      xd1 = rbind(kronecker(jp,(diag(sig)/lamda)), 
                  matrix(0, nrow = n, ncol = (n*p)))
    }
  }
  
  # Get additional dummy matrices - equation (9) Banbura et al. 2007
  if (tau >0){
    if (epsilon > 0){
      delta = c(delta)
      mu = c(mu)
      yd2 = diag(delta*mu)/tau
      xd2 = cbind(kronecker(matrix(1,nrow=1,ncol=p),yd2),0)
    }else{
      yd2 = diag(delta*mu)/tau
      xd2 = kronecker(matrix(1,nrow=1,ncol=p),yd2)
    }
    
  }
  y = rbind(yd1, yd2)
  x = rbind(xd1, xd2)
  return(list(y = y, x = x))
}

#create_dummiestrendx###########################################################
create_dummiestrendx <- function(lamda,tau,delta,epsilon,p,mu,sig,n){
  # Creates matrices of dummy observations, following Banbura et al. 2007
  # lamda tightness parameter
  # tau  prior on sum of coefficients
  # delta prior mean for VAR coefficients
  # epsilon tigtness of the prior around constant
  # mu sample mean of the data
  # sig vector of AR residual variances for the data, input from getdummies()
  
  #Initialize vectors
  yd1 = c()
  yd2 = c()
  xd1 = c()
  xd2 = c()
  
  #Get dummy matrices in equation (5) of Banbura et al. 2007
  if (lamda > 0){
    if (epsilon > 0){
      #yd
      yd1 = rbind(diag(sig*delta)/lamda, 
                  matrix(0,nrow =(n*(p-1)), ncol = n))
      yd1 = rbind(yd1, diag(sig))
      yd1 = rbind(yd1, 0, 0)
      
      #xd
      jp = diag(1:p)
      xd1 = cbind(kronecker(jp,(diag(sig)/lamda)), 
                  matrix(0, nrow = (n*p), ncol = 2))
      xd1 = rbind(xd1, matrix(0, nrow = n, ncol = (n*p)+2))
      xd1 = rbind(xd1, cbind(matrix(0, nrow = 1, ncol = (n*p)), epsilon[1], 0))
      xd1 = rbind(xd1, cbind(matrix(0, nrow = 1, ncol = (n*p)), 0 , epsilon[1]))
    }else{
      #yd
      yd1 = rbind(diag(sig*delta)/lamda, 
                  matrix(0,nrow =(n*(p-1)), ncol = n))
      yd1 = rbind(yd1, diag(sig))
      
      #xd
      jp = diag(1:p)
      xd1 = rbind(kronecker(jp,(diag(sig)/lamda)), 
                  matrix(0, nrow = n, ncol = (n*p)))
      rm(jp)
    }
  }
  
  # Get additional dummy matrices - equation (9) Banbura et al. 2007
  if (tau >0){
    deltav = as.vector(delta)
    muv = as.vector(mu)
    if (epsilon >0){
      yd2 = diag(deltav*muv)/tau
      xd2 = cbind(kronecker(matrix(1,nrow=1,ncol=p),yd2),
                  matrix(0, nrow=n, ncol=2))
      
    }else{
      yd2 = diag(deltav*muv)/tau
      xd2 = kronecker(matrix(1,nrow=1,ncol=p), yd2)
    }
  }
  rm(list = c("deltav", "muv"))
  y = rbind(yd1, yd2)
  x = rbind(xd1, xd2)
  
  rm(list = c("yd1", "xd1", "yd2", "xd2"))
  return(list(y = y, x = x))
}

#extract########################################################################
extract <- function(data, k){
  # This function extracts the first k principal components from a (txn) matrix
  # 'data' and returns the factors (fac, (txk)) as well as the normalised 
  # 'loading  matrix (lam, (nxk)) - i.e, the matrix composed of eigenvectors
  # of the data's covariance matrix such that the principal components are equal
  # to: fac= = data*lam
  #
  # Note 1: we normalise to ensure that t(lam)*lam/n = I
  # Note 2: Note: The function does not check whether the means of the variables
  # are zero, if they are not it will return incorrect results
  
  #matrix dimensions
  t = nrow(data); n = ncol(data)
  
  # Assuming the empirical mean of the columns are 0, the empirical covariance 
  # matrix is:
  xx = t(data)%*%data
  eval = eigen(xx)$values; evec = eigen(xx)$vectors
  
  # Obs: In the original function, they reorder the eigenvalues and eigenvectors
  # in the descending order. Here we do not need to do such step, since the 
  # eigen() function already return both values and vectors in the descending 
  # order. See function documentation
  
  # take the first k eigenvectors and normalise by square root of data series
  lam = sqrt(n)*evec[ ,1:k]
  
  # Compute the principal component vectors from standard formula. Note here
  # all the principal components are 'shortned' by the division
  fac = data%*%lam/n
  
  return(list(fac = fac, lam = lam))
}
#get_dummies####################################################################
get_dummies <- function(LAMDAPEVIEWS,TAUPEVIEWS,EPSILONPEVIEWS,Y,L,EX,LH,EPSILONH,RW){
  # Initialize the AR process to get the minessota prior with dummy variables
  # following Banbura et al. 2007
  # !!! Requires pracma library for Moore Pseudo-inverse pinv()
  
  library("pracma")
  
  #N = ncol(Y)
  #lamdaP=LAMDAPEVIEWS       # tightness of the priors on the first lag
  #tauP=TAUPEVIEWS           # tightness of the priors on sum of coefficients
  #epsilonP=EPSILONPEVIEWS   # tightness of the prior on the constant
  #epsilonH=EPSILONH 
  #muP=t(colMeans(Y))
  sigmaP= numeric(N)
  deltaP= numeric(N)
  
  # AR step
  for (i in 1:ncol(Y)){
    ytemp = matrix(Y[,i], nrow=nrow(Y), ncol=1)
    ytemplag = matrix(c(NA, ytemp[(1:(nrow(ytemp) - 1))]), 
                      nrow = nrow(Y), ncol = 1)
    xtemp = cbind(ytemplag, 1)
    ytemp = ytemp[2:nrow(ytemp),]
    xtemp = xtemp[2:nrow(xtemp),]
    btemp = solve(t(xtemp)%*%xtemp)%*%t(xtemp)%*%ytemp
    etemp = ytemp -xtemp%*%btemp
    rm(xtemp)
    stemp = t(etemp)%*%etemp/length(ytemp)
    rm(list = c("ytemp", "etemp"))
    if (btemp[1] > 1){
      btemp[1] = 1
    }
    deltaP[i] = btemp[1]
    if (RW){
      deltaP[i] = 1
    }
    rm(btemp)
    sigmaP[i] = sqrt(stemp + 0i) #+0i so sqrt wont return Nans for negative values
    rm(stemp)
  }
  
  dummies = create_dummies(LAMDAPEVIEWS, TAUPEVIEWS, deltaP, EPSILONPEVIEWS, L, 
                           t(colMeans(Y)), sigmaP, ncol(Y), LH, EPSILONH)
  xd = dummies$x
  yd = dummies$y
  rm(dummies)
  b0 = solve(t(xd)%*%xd)%*%t(xd)%*%yd
  
  #b01= matrix(b0,nrow=(N*L + EX), ncol=N)
  e0 =yd - xd%*%b0
  S = diag(sigmaP^2)
  s0 = kronecker(S,pinv(t(xd)%*%xd))
  rm(S)
  return(list(yd = yd, xd = xd, b0=b0, s0=s0))
}

#getinitialvol##################################################################
getinitialvol <- function(data0,REPS,BURN,T0,L){
  library("pracma")
  
  # A gibbs sampling procedure to get priors for volatility 
  N = ncol(data0)
  Y = data0
  trend = as.matrix((0:(nrow(Y) - 1)), ncol=1)
  #take lags
  X = c()
  for (j in 1:12){                 
    lags = lag0(Y, k = j)
    X = cbind(X, lags)
  }
  
  X = cbind(X,1,trend)
  
  Y = Y[(L+1):nrow(Y),]
  X = X[(L+1):nrow(X),]
  T = nrow(X)
  
  # initial training sample
  Y0 = Y[1:T0,]
  X0 = X[1:T0,]
  # Prior for VAR coefficients
  lamdaP = 0.1 #Tightness of the priors on the first lag
  tauP = 10*lamdaP #tightness of the priors on sum of coefficients
  epsilonP = 1/1000 #tightness of the prior on the constante
  muP=t(colMeans(Y0))
  sigmaP= numeric(N)
  deltaP= numeric(N)
  e0 = c()
  for (i in 1:N){
    ytemp = matrix(Y0[,i], nrow=nrow(Y), ncol=1)
    ytemplag = matrix(c(NA, ytemp[(1:(nrow(ytemp) - 1))]), 
                      nrow = nrow(Y), ncol = 1)
    xtemp = cbind(ytemplag, 1)
    ytemp = ytemp[2:nrow(ytemp),]
    xtemp = xtemp[2:nrow(xtemp),]
    #btemp = solve(t(xtemp)%*%xtemp)%*%t(xtemp)%*%ytemp
    btemp = mldivide(xtemp,ytemp)
    etemp = ytemp -xtemp%*%btemp
    stemp = t(etemp)%*%etemp/length(ytemp)
    if (btemp[1] > 1){
      btemp[1] = 1
    }
    deltaP[i] = btemp[1]
    sigmaP[i] = sqrt(stemp) #+0i so sqrt wont return Nans for negative values
    e0 = cbind(e0, etemp)
  }
  
  # Dummy data to implement priors
  dummiestrendx = create_dummiestrendx(lamdaP,tauP,deltaP,epsilonP,L,muP,sigmaP,N)
  xd = dummiestrendx$x
  yd = dummiestrendx$y
  
  # Prior mean and variance using the dummy data
  # B0 = solve(t(xd)%*%xd)%*%t(xd)%*%yd
  B0 = mldivide(xd,yd)
  E0 = yd - xd%*%B0
  SIGMA0 = kronecker((t(E0)%*%E0),pinv(t(xd)%*%xd)) #prior variance
  B0 = matrix(B0, nrow = length(B0), ncol = 1)
  # prior for the elements of the A matrix
  s0 = (t(e0)%*%e0)/T0
  C0 = chol(s0)
  C0 = t(solve(C0/rep(diag(C0),N)))
  SC0 = 10 #prior variance
  
  # prior for the initial condition of the stochastic volatility
  MU0 = log(diag(s0)) # prior mean
  SV0 = 1             # prior variance
  
  # remove training sample
  Y = Y[(T0+1):nrow(Y),]
  X = X[(T0+1):nrow(X),]
  T = nrow(Y)
  
  #rough guess for stochastic volatility
  hlast = (diff(Y)^2) + 0.0001
  hlast = rbind(hlast[1:2,],hlast) #rough initial guess for svol
  errors = diff(Y) # initial guess for VAR residuals
  errors = rbind(as.vector(errors[1,]),errors[(1:nrow(errors)),])
  g = matrix(1,nrow = N, ncol=1) #rough guess for the variance of the trans. eq.
  g0 =  0.01^2 # prior scale paramater for inverse gamma
  Tg0 = 1 # prior degrees of freedom
  # SV=zeros(REPS-BURN,T+1,N) -> Mumtaz didnt use this variable
  h01 = hlast
  h0 = hlast
  errors0 = errors
  for (i in 1:50){
    for (j in 1:N){
      h01[,j] = getsvol(h0[,j],g[j],log(s0[j,j]),10,errors[,j])
    }
    h0 = h01
  }
  
  hlast = h01
  B0h = rbind(0.7,0)
  Sigma0h = diag(c(0.5,0.00005))
  MUSVOL = numeric(N)
  FSVOL = numeric(N)
  
  beta20 = as.vector(mldivide(X,Y))
  
  # pre allocate variables to store values
  out1 = c()
  out2 = array(0, dim = c((REPS - BURN), T+1, N))
  out3 = array(0, dim = c((REPS - BURN), T, (N*(N-1))/2))
  out4 = array(0, dim = c((REPS - BURN), N))
  out5 = array(0, dim = c((REPS - BURN), N*2))
  
  jgibbs = 1
  for (igibbs in 1:REPS){
    print(paste("Gibbs iteration number: ", igibbs))
    #Step 1 of the Gibbs Algorithm: sample the A matrix
    amatx = c()
    
    for (j in 2:N){
      # v2=-a1*v1+sqrt(exp(h2))*e2
      ytemp = errors[,j]
      xtemp = errors[,(1:(j-1))]*(-1)
      ytemp = ytemp/(sqrt(hlast[2:nrow(hlast),j])) # remove heteroscedasticity
      xtemp = xtemp/(rep(sqrt(hlast[2:nrow(hlast),j]), n=(j-1))) # remove heteroscedasticity
      
      # prior means and variance
      A0 = t(t(C0[j, 1:(j-1)]))
      
      if ((j-1) == 1){
        V00 = diag(as.matrix(abs(A0)))*SC0
      }else{
        V00 = diag(c(abs(A0)))*SC0
      }
      
      atemp = getregx(ytemp,xtemp,A0,V00,1)
      
      amatx = cbind(amatx, t(atemp))
    }
    
    #Step 2 of the algorithm: Sample stochastic volatilities
    A = chofac(N, t(amatx))
    epsilon = errors%*%t(A) # orthogonal residuals
    
    # sample stochastic vol for each epsilon using the MH algorithm in Jaquier, 
    # Polson and Rossi
    hnew = c()
    for (i in 1:N){
      htemp= getsvolx(hlast[,i],g[i],MU0[i],SV0,epsilon[,i],MUSVOL[i],FSVOL[i])
      hnew = cbind(hnew,htemp)
    }
    hlast=hnew
    #draw AR parameters
    gerrors = c()
    for (jj in 1:N){
      
      ytemp = log(hlast[,jj])
      xtemp = cbind(matrix(c(0,diff(ytemp)), nrow =  length(ytemp))
                    , matrix(1, nrow = length(ytemp), ncol=1))
      ytemp = ytemp[2:nrow(xtemp)]
      xtemp = xtemp[2:nrow(xtemp),]
      
      MM = solve(solve(Sigma0h) + (1/g[jj])*(t(xtemp)%*%xtemp))%*%(solve(Sigma0h)%*%B0h + (1/g[jj])*t(xtemp)%*%ytemp)
      VV = solve(solve(Sigma0h) + (1/g[jj])*(t(xtemp)%*%xtemp))
      
      chck=-1
      while(chck<0){
        BB = MM + t(matrix(rnorm(2), ncol=2)%*%chol(VV))
        ee = max(abs(BB[1]))
        if (ee <=1){
          chck=1
        }
      }
      FSVOL[jj] = BB[1]
      MUSVOL[jj] = BB[2]
      gerrors = cbind(gerrors, (ytemp - xtemp%*%BB))
    }
    
    # Step 3 of the algorithm: Sample g from the IG distribution
    for (i in 1:N){
      g[i] = IG(Tg0,g0,gerrors[,i])
    }
    carterkohn = carterkohnvar(Y,X,0,invpd(A),hlast,t(B0),SIGMA0,L,1,50,2)
    problem = carterkohn$problem
    beta2 = carterkohn$beta2
    if (problem){
      beta2 = beta20
    }else{
      beta20 = carterkohn$beta2
    }
    
    if (igibbs > BURN){
      out1 = rbind(out1, beta2)
      out2[jgibbs, 1:(T+1), 1:N] = hlast[1:nrow(hlast),]
      out3[jgibbs, 1:T, 1:((N*(N-1))/2)]= repmat(amatx, T, 1)
      out4[jgibbs, 1:N] = t(g)
      out5[jgibbs, 1:(N*2)] = cbind(t(FSVOL),t(MUSVOL))
      jgibbs = jgibbs + 1
    }
  }
  
  outh = drop(rowMeans(out2))
  outf = colMeans(out5)
  outg = colMeans(out4)
  tmp = drop(rowMeans(out3))
  outA = tail(tmp, n=1)
  
  
  return(list(outh = outh, outf = outf, outg = outg, outA = outA))
}

#getqr##########################################################################
getqr <- function(a){
  # Returns a modifiedQR decomposition of a matrix a, where diagonal elements of
  # the R matrix are all positive
  temp = qr(a)
  q = qr.Q(temp)
  r = qr.R(temp)
  for (i in nrow(q)){
    if (r[i,i] < 0) {
      q[,i] = -1*q[,i]
    }
  }
  return(q)
}
#getregx#########################################################################
getregx <- function(Y,X,B0,SIGMA0,sigma2){
  isigma0 = invpd(SIGMA0)
  isigmatemp = isigma0 
  # checks if isigma0 is a vector of size 1, if yes change it to vector so R
  # does not freak out
  if(length(SIGMA0) == 1){
    isigmatemp = c(invpd(SIGMA0))
  }
  isigma2 = 1/sigma2
  xx = t(X)%*%X #the opposite of what Mumtaz does, but in this produce the actual
  #result in the form needed
  V = invpd(isigmatemp + isigma2*xx)
  M = V%*%(isigma0%*%B0 + isigma2*t(X)%*%Y)
  rm(list = c("isigma0", "isigma2"))
  
  if(typeof(ncol(X)) == "NULL"){
    bdraw = M + t(matrix(rnorm(1), ncol = 1)%*%chol(V))
  }else{
    bdraw = M + t(matrix(rnorm(ncol(X)), ncol = ncol(X))%*%chol(V))
  }
  
  return(bdraw)
}

#getsvol########################################################################
getsvol <- function(hlast,g,mubar,sigmabar,errors){
  
  T = length(errors)
  hnew = numeric(T)
  
  i=0
  # Time period 0 ###
  
  hlead = hlast[i+1]
  ss = sigmabar*g/(g + sigmabar) # variance
  mu = ss*((mubar/sigmabar) + (log(hlead)/g)) # mean
  #draw from lognormal using mu and ss
  h = exp(mu + ((ss^0.5)*rnorm(1)))
  hnew[i + 1] = h
  
  # time period 1 to t-1 ###
  for (i in 2:T){
    hlead = hlast[i+1]
    hlag = hnew[i-1]
    yt = errors[i-1]
    
    # mean and variance of the proposal log normal density
    mu = (log(hlead) + log(hlag))/2
    ss = g/2
    
    # candidate draw from lognormal
    htrial = exp(mu + (ss^0.5)*rnorm(1))
    
    # acceptance probability in logs
    lp1 = -0.5*log(htrial) - (yt^2)/(2*htrial) # numerator
    lp0 = -0.5*log(hlast[i]) - (yt^2)/(2*hlast[i]) # denominator
    accept = min(1, exp(lp1 - lp0)) # ensure accept <=1
    
    u = rand(1,1)
    if (u <= accept){
      h = htrial
    }
    else{
      h = hlast[i]
    }
    hnew[i] = h
  }
  
  # Time period T ###
  i = T+1
  yt = errors[i-1]
  hlag = hnew[i-1]
  
  #mean and variance of the proposal density
  mu = log(hlag)
  ss = g
  
  # candidate draw from lognormal
  htrial = exp(mu + ((ss^0.5)*rnorm(1)))
  
  # acceptance probability
  lp1 = -0.5*log(htrial) - (yt^2)/(2*htrial) # numerator
  lp0 = -0.5*log(hlast[i]) - (yt^2)/(2*hlast[i]) # denominator
  accept = min(1, exp(lp1 - lp0)) # ensure accept <=1
  
  u = rand(1,1)
  if (u <= accept){
    h = htrial
  }
  else{
    h = hlast[i]
  }
  hnew[i] = h
  
  rm(list = c("h", "T", "htrial", "u", "accept", "lp0", "lp1", "ss", "mu", "yt"))
  return(hnew)
}

#getsvolx#######################################################################
getsvolx <- function(hlast,g,mubar,sigmabar,errors,alpha,delta){
  T = length(errors)
  hnew = numeric(T)
  
  i = 0
  # Time period 0
  
  hlead = hlast[i+1]
  ss = sigmabar*g/(g + ((delta^2)*sigmabar)) # variance
  mu = ss*((mubar/sigmabar) + delta*(log(hlead)-alpha)/g)
  
  #draw from lognormal using mu and ss
  h = exp(mu + ((ss^0.5)*rnorm(1)))
  hnew[i+1] = h
  
  # time period 1 to T-1
  for (i in 2:T){
    hlead = hlast[i+1]
    hlag = hnew[i-1]
    yt = errors[i-1] 
    
    #mean and variance of the proposal log normal density
    mu = alpha*(1-delta) + delta*(log(hlead) + log(hlag))/(1 + delta^2)
    ss = g/(1 + delta^2)
    
    #candidate draw from lognormal
    htrial = exp(mu + (ss^0.5)*rnorm(1))
    
    #acceptance probability in logs
    lp1 = -0.5*log(htrial) - (yt^2)/(2*htrial) #numerator
    lp0 = -0.5*log(hlast[i]) - (yt^2)/(2*hlast[i]) #denominator
    accept = min(1, exp(lp1 - lp0)) #ensure acceptance <=1
    
    u = rand(1,1)
    if (u <= accept){
      h = htrial
    }
    else{
      h = hlast[i]
    }
    hnew[i]= h
  }
  # time period T
  i = T + 1
  yt = errors[i-1]
  hlag = hnew[i-1]
  
  #mean and variance of the proposal density
  mu = alpha + delta*log(hlag)
  ss = g
  
  #candidate draw from lognormal
  htrial = exp(mu + (ss^0.5)*rnorm(1))
  
  # acceptance probability
  lp1 = -0.5*log(htrial) - (yt^2)/(2*htrial)
  lp0 = -0.5*log(hlast[i]) - (yt^2)/(2*hlast[i])
  accept = min(1, exp(lp1 - lp0)) # ensure acceptance<=1
  
  u = rand(1,1)
  if (u <= accept){
    h = htrial
  }
  else{
    h = hlast[i]
  }
  hnew[i] = h
  
  rm(list = c("h", "T", "htrial", "u", "accept", "lp0", "lp1", "ss", "mu", "yt"))
  return(hnew)
}

#invpd##########################################################################
invpd <- function(x){
  library("pracma")
  as.matrix(x)
  if (typeof(ncol(x)) == "NULL"){
    temp = 1
  }else{
    temp = diag(1, nrow = ncol(x))
  }
  
  out = mldivide(x,temp)
  
  rm(temp)
  return(out)
}

#IG#############################################################################
IG <- function(v0,d0,x){
  # return posterior draw - v - from a inverse gamma with prior degrees of 
  # freedom v0/2 and scale parameter d0/2. The posterior values are v1 and d1,
  # respectively. x is a vector of innovations
  
  # The simulation method follows bauwens et al. p 317 IG2(s,v):
  #   simulate x =chisquare(V)
  #   deliver = s/x
  T = length(x)
  v1 = v0 + T
  d1 = d0 + (t(x)%*%x)
  z = rnorm(v1)
  x = t(z)%*%z
  v = d1/x
  
  rm(list = c("T", "v1", "d1", "z", "x"))
  return(v)
}
#IWPQ###########################################################################
IWPQ <- function(v, ixpx){
  # returns inverse wishart prior
  k = nrow(ixpx)
  z = matrix(0, nrow = v, ncol = k)
  mu = matrix(0, nrow = k, ncol = 1)
  for (i in 1:v){
    z[i,] = t(t(chol(ixpx))%*%rnorm(k))
  }
  return(solve(t(z)%*%z))
}
#irfsim#########################################################################
irfsim <- function(b,n,l,v,s,t){
  # b = VAR coefficients
  # n = number of variables
  # l = lag length
  # v = A0 matrix
  # s =  shock vector
  # t = horizon
  
  e = matrix(0, nrow = (t+1), ncol = n)
  e[(l+1),] = s
  y = matrix(0, nrow = (t+1), ncol = n)
  
  for (k in (l+1):t){
    x = c()
    for (i in 1:l){
      for (j in 1:n){
        x = cbind(x, y[(k-i),j])
      }
    }
    y[k,] = (cbind(x,0)%*%b) + (e[k,]%*%v)
  }
  y = y[(l+1):(nrow(y) - 1),]
  return(y)
}
#lag0###########################################################################
lag0 <- function(x, k) {
  # ensure 'x' is a matrix
  stopifnot(is.matrix(x))
  if (k == 0)
    return(x)
  na <- matrix(0, nrow=abs(k), ncol=ncol(x))
  if (k > 0) {
    nr <- nrow(x)
    # prepend NA and remove rows from end
    return(rbind(na, as.matrix(x[1:(nr - k),])))
  } else {
    # append NA and remove rows from beginning
    return(rbind(x[-1:k,], na))
  }
}
#stability######################################################################
stability <- function(beta,n,l,ex){
  # coef -> (n*1 + 1)Xn matrix with the coef from the VAR
  # 1 -> number of lags
  # n -> number of endogenous variables
  # FF -> matrix with all coef
  # returns S -> a dummy var. If S = 1 (TRUE), this means instability
  
  FF = matrix(0, nrow =(n*l), ncol = (n*l))
  FF[(n+1):(n*l),1:(n*(l-1))] = diag(1, nrow = (n*(l-1)), ncol = (n*(l-1)))
  
  temp = matrix(beta, nrow = ((n*l) + ex), ncol = n)
  temp = t(temp[1:(n*l), 1:n])
  FF[1:n,1:(n*l)] = temp
  ee = max(abs(eigen(FF)$values))
  S = (ee > 1)
  
  rm(list = c("FF", "ee", "temp"))
  return(S)
}
#stabilityABE###################################################################
stabilityABE <- function(beta,n,l){
  # This the function he uses in Applied Bayesian Econometric Book
  # 1 -> number of lags
  # n -> number of endogenous variables
  # FF -> matrix with all coef
  # returns S -> a dummy var. If S = 1 (TRUE), this means instability
  
  FF = matrix(0, nrow =(n*l), ncol = (n*l))
  FF[(n+1):(n*l),1:(n*(l-1))] = diag(1, nrow = (n*(l-1)), ncol = (n*(l-1)))
  
  temp = matrix(beta, nrow = ((n*l) + 1), ncol = n)
  temp = t(temp[2:((n*l) + 1), 1:n])
  FF[1:n,1:(n*l)] = temp
  ee = max(abs(eigen(FF)$values))
  S = (ee > 1)
  
  rm(list = c("FF", "ee", "temp"))
  return(S)
}
#standardise####################################################################
standardise <- function(x){
  # takes a matrix x and return it in standardise
  t = nrow(x); n = ncol(x)
  
  # matrice m and s, with the same dimension as x, with the mean and sd of 
  # every column of x
  sd_x = apply(x, 2, sd)
  mean = colMeans(x)
  m = matrix(rep(mean, t), nrow = t, byrow = TRUE)
  s = matrix(rep(sd_x, t), nrow = t, byrow = TRUE)
  
  # Subtract the mean and divide by the standard deviation
  y=(x-m)/s
  
  return(y)
}