clear
addpath('functions'); %this line adds functions to take lags etc
%an AR 2 model for US inflation
%load inflation data
Y=xlsread('\data\inflation.xls');
T=rows(Y);
X=[ones(T,1) lag0(Y,1) lag0(Y,2)];

%remove missing obs
Y=Y(3:end);
X=X(3:end,:);
T=rows(X);
%step 1 set priors and starting values
%priors for B
B0=[0;0;0];
Sigma0=eye(3);
%priors for sigma2
T0=1;
D0=0.1;
%starting values
B=B0;
sigma2=1;

reps=5000;   %total numbers of Gibbs iterations
burn=4000;   %percent of burn-in iterations
out1=[];
out2=[];
out3=[];
for i=1:reps
%step 2 Sample B conditional on sigma N(M*,V*)
M=inv(inv(Sigma0)+(1/sigma2)*(X'*X))*(inv(Sigma0)*B0+(1/sigma2)*X'*Y); 
V=inv(inv(Sigma0)+(1/sigma2)*(X'*X));
chck=-1;
while chck<0                     %check for stability
B=M+(randn(1,3)*chol(V))';
b=[B(2) B(3);1   0];
ee=max(abs(eig(b)));
if ee<=1
    chck=1;
end
end
%step 3 sample sigma2 conditional on B from IG(T1,D1);
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
    
    %compute forecast for 2 years
    yhat=zeros(14,1);
    yhat(1:2)=Y(end-1:end); %starting values
    cfactor=sqrt(sigma2);  %standard deviation of the shocks
    for m=3:14
        yhat(m)=[1 yhat(m-1) yhat(m-2)]*B+(randn(1,1)*cfactor);
    end
    %save
    out3=[out3 [Y;yhat(3:end)]];    
end
end

%plot marginal posterior distributions
figure(1)
subplot(2,2,1);
hist(out1(:,1),50);
axis tight
title('Constant');


subplot(2,2,2);
hist(out1(:,2),50);
axis tight
title('AR(1) coefficient');

subplot(2,2,3);
hist(out1(:,3),50);
axis tight
title('AR(2) coefficient');


subplot(2,2,4);
hist(out2(:,1),50);
axis tight
title('\sigma^{2}');


%compute mean of the marginal posterior distribution of B
MB=mean(out1);
%compute standard error
VB=std(out1);
%compute 95% error band
EB=prctile(out1,[5 95]);

%plot forecast distribution
figure(2)
TT=1947.25:0.25:2012.5;
outx=prctile(out3',[10 20 30 40 50 60 70 80 90])'; %here prctile calculates percentiles of the forecast distribution
plot(TT,outx);

xlim([2000 2013])


