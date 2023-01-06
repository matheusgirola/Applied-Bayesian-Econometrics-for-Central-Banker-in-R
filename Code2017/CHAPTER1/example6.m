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

%analytical computation of the marginal likelihood
mlika=mlikols(B0,Sigma0,T0,D0,Y,X);
disp('Analytical log Marginal Likelihood');
disp(log(mlika));

sigma2=1;
reps=15000;   %total numbers of Gibbs iterations
burn=4000;   %percent of burn-in iterations
out1=[];
out2=[];
for i=1:reps

%Sample B conditional on sigma N(M*,V*)
M=inv(inv(Sigma0)+(1/sigma2)*(X'*X))*(inv(Sigma0)*B0+(1/sigma2)*X'*Y); 
V=inv(inv(Sigma0)+(1/sigma2)*(X'*X));

B=M+(randn(1,2)*chol(V))';

%sample sigma2 conditional on B from IG(T1,D1);
%compute residuals
resids=Y-X*B;
%compute posterior df and scale matrix
T1=T0+T;
D1=D0+resids'*resids;
%draw from IG
z0=randn(T1,1);
z0z0=z0'*z0;
sigma2=D1/z0z0;

if i>burn
    out1=[out1;B'];
    out2=[out2;sigma2];
end
end

%calculate the marginal likelihood using Chib's method
%posterior mean
bstar=mean(out1);
sigmastar=mean(out2);
Hstar=mean(1./out2);
% Step1 evaluate the prior distributions at the posterior mean
%P(B)~N(B0,Sigma0) in logs
Pb=log(mvnpdf(bstar',B0,Sigma0));
%P(1/sigma2)~Gamma(D0,T0) in logs
PH=gampdf1(T0,D0,Hstar);
%Step 2 evaluate the log likelihood
loglik=-(T/2)*log(2*pi*sigmastar)-0.5*(((Y-X*bstar')'*(Y-X*bstar'))/sigmastar);

%step 3 evaluate the posterior density
%H(bstar,sigmastar)=H(bstar\sigmastar)*H(sigmastar)
%step 3a H(bstar\sigmastar)~N(M1,V1) in logs
M1=inv(inv(Sigma0)+(1/sigmastar)*(X'*X))*(inv(Sigma0)*B0+(1/sigmastar)*X'*Y); 
V1=inv(inv(Sigma0)+(1/sigmastar)*(X'*X));
H1=log(mvnpdf(bstar',M1,V1));

%step 3b evaluate H(sigmastar) using the Gibbs draws
H2log=[];
for i=1:size(out1,1)
    bgibbs=out1(i,:)';
    res=Y-X*bgibbs;
   H2i=gampdf1(rows(res)+T0,(res'*res)+D0,Hstar);
    H2log=[H2log;H2i];
end

%take exponential and mean in a way that prevents underflow
factor=max(H2log);
H2exp=exp(H2log-factor);
H2mean=log(mean(H2exp))+factor;

%calculate marginal likelihood
mlik=loglik+Pb+PH-H1-H2mean;

disp('Chib log Marginal Likelihood');
disp(mlik);


