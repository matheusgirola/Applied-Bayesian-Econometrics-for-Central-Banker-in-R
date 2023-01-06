clear;
clc
addpath('functions')
%generate artificial data
T=100;
X=[ones(T,1) randn(T,1)];
btrue=[1;0.5];
sigmatrue=0.2;
Y=X*btrue+randn(T,1)*sqrt(sigmatrue);

%set priors
T0=3;
D0=2.5;

B0=zeros(2,1);
Sigma0=eye(2)*(4);

%step 2 set SIGMA matrix via OLS estimation
yols=Y;
xols=X;
bols=inv(xols'*xols)*(xols'*yols);
eols=yols-xols*bols;
sols=((eols'*eols)/T);
vols=sols*inv(xols'*xols);
K=0.1;
P=eye(3);                    %this is the variance of the metropolis hastings random walk based partly on OLS estimates
P(1,1)=(vols(1,1));
P(2,2)=(vols(2,2));
P(3,3)=0.1;


%analytical computation of the marginal likelihood
mlika=mlikols(B0,Sigma0,T0,D0,Y,X);
disp('Analytical log Marginal Likelihood');
disp(log(mlika));

sigma2=1;
Gammaold=[0;0;1]; %starting values
%compute posterior
posteriorOLD=postols(Y,X,Gammaold,B0,Sigma0,T0,D0);
reps=15000;   %total numbers of MH iterations
burn=4000;   %percent of burn-in iterations
outpost=[];  %will hold posterior
outparam=[]; %will hold parameters
naccept=0;
for i=1:reps

Gammanew=Gammaold+(randn(1,3)*chol(P*K))';
  sigma2=Gammanew(3);
    if sigma2<0
        posteriorNEW=-1000000;
    else
        posteriorNEW=postols(Y,X,Gammanew,B0,Sigma0,T0,D0);
    end
        accept=min([exp(posteriorNEW-posteriorOLD);1]);   %min(accept,1)
    
    u=rand(1,1);  %random number from the uniform dist
    
    if u<accept
        Gammaold=Gammanew;  %accept draw
        naccept=naccept+1;  %count number of acceptances
        posteriorOLD=posteriorNEW;
      
    end

if i>burn
    outpost=[outpost;posteriorOLD];
     outparam=[outparam;Gammaold'];
end
end

%calculate the marginal likelihood using Gelfand and Dey method
%posterior mean and variance
pmean=mean(outparam);
pvar=cov(outparam);
lpost_mode=max(outpost);
p=0.1; %critical value of the Chi-squared distribution
npara=size(outparam,2); %number of parameters
critval = chi2inv(p,npara);
    tmp = 0;
    for i = 1:size(outparam,1);
        deviation  = (outparam(i,:)-pmean)*inv(pvar)*((outparam(i,:)-pmean))';
        if deviation <= critval;
            lftheta = -log(p)-(npara*log(2*pi)+log(det(pvar))+deviation)/2;
            tmp = tmp + exp(lftheta - outpost(i)+lpost_mode);
        end;    
    end;
mlik=lpost_mode-log(tmp/size(outparam,1));

disp('Gelfand and Dey log Marginal Likelihood');
disp(mlik);


