clear
clc
addpath('functions');
REPS=15000;
BURN=10000;
[data,names]=xlsread('\data\usdata1.xls'); %load US data
N=cols(data);

L=2; %lag length of the VAR
Y=data;
%take lags
X=[];
for j=1:L
X=[X lag0(data,j) ];
end
X=[X ones(rows(X),1)];
Y=Y(L+1:end,:);
X=X(L+1:end,:);
T=rows(X);
%priors for VAR coefficients
lamdaP=1;  %This controls the tightness of the priors on the first lag
tauP=0;  % this controls the tightness of the priors on sum of coefficients 0 means not applied
epsilonP=1;  % this controls tightness of the prior on the constant
muP=mean(Y)';
sigmaP=[];
deltaP=[];
for i=1:N
    ytemp=Y(:,i);
    xtemp=[lag0(ytemp,1) ones(rows(ytemp),1)];
    ytemp=ytemp(2:end,:);
    xtemp=xtemp(2:end,:);
    btemp=xtemp\ytemp;
    etemp=ytemp-xtemp*btemp;
    stemp=etemp'*etemp/rows(ytemp);
    deltaP=[deltaP;btemp(1)];
    sigmaP=[sigmaP;stemp];
end
%dummy data to implement priors see http://ideas.repec.org/p/ecb/ecbwps/20080966.html
[yd,xd] = create_dummies(lamdaP,tauP,deltaP,epsilonP,L,muP,sigmaP,N);

%yd and xd are the dummy data. Append this to actual data
  Y0=[Y;yd];
  X0=[X;xd];
  
 %compute the marginal likelihood analytically for comparison
 temp1=mlikvar1(Y,X,yd,xd);
  %conditional mean of the VAR coefficients
  mstar=vec(X0\Y0);  %ols on the appended data
  xx=X0'*X0;
  ixx=xx\eye(cols(xx));  %inv(X0'X0) to be used later in the Gibbs sampling algorithm
  sigma=eye(N); %starting value for sigma
  out1=zeros(REPS-BURN,N*(N*L+1),1);
  out2=zeros(REPS-BURN,N,N);
  jj=1;
  for i=1:REPS
       vstar=kron(sigma,ixx);
       beta=mstar+(randn(1,N*(N*L+1))*chol(vstar))';
       
       %draw covariance
       e=Y0-X0*reshape(beta,N*L+1,N);
    scale=e'*e;
    sigma=IWPQ(T+rows(yd),inv(scale));
    
    if i>=BURN
        out1(jj,:,:)=beta;
        out2(jj,:,:)=sigma;
        jj=jj+1;
    end
  end
  
  betam=squeeze(mean(out1,1));
  sigmam=squeeze(mean(out2,1));
  
  %evaluate priors
  b0=vec(xd\yd);
  b01=reshape(b0,N*L+1,N);
  e0=yd-xd*b01;
  S=e0'*e0;
  
  %evaluate log prior distribution for VAR coefficients
  bp=multivariatenormal(betam',b0,kron(S,pinv(xd'*xd)));
  %evaluate log prior for VAR covariance
  sp= invwishpdf(sigmam,S,size(yd,1)-size(xd,2));
  %evaluate log likelihood
  lik=loglik(reshape(betam,N*L+1,N),sigmam,Y,X);
  %evaluate H(Bstar\sigmastar);
  vstar1=kron(sigmam,ixx);
  H1=multivariatenormal(betam',mstar,vstar1);
  %evaluate H(sigmastar\beta[j])
  H2i=[];
  for j=1:size(out1,1)
      betaj=out1(j,:);
      e=Y0-X0*reshape(betaj,N*L+1,N);
      scale=e'*e;
      H2i= [H2i;invwishpdf(sigmam,scale,size(Y0,1))];
  end
  %take mean taking care of possible underflow/overflow with exp
  factor=max(H2i);
  H2=exp(H2i-factor);
  H2m=mean(H2);
  H2m=log(H2m)+factor;
  
  %marginal lik
  mlik=lik+bp+sp-H1-H2m;
  
  disp('Analytical log Marginal Likelihood')
  disp(temp1);
  
  disp('Chib log Marginal Likelihood')
  disp(mlik);
  
  
  
