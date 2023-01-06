clear
addpath('functions');
%an AR 2 model for US inflation with autocorrelated AR(1) disturbances 
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
%priors for rho
rho0=0;
Sigma0r=1;

%starting values
B=B0;
rho=rho0;
sigma2=1;

reps=25000;
burn=24000;
out1=[]; %will save the inflation forecast
out2=[];
out3=[];
out4=[];
for i=1:reps

%step 2 Sample B conditional on sigma N(M*,V*)
%remove serial correlation
ystar=Y-lag0(Y,1)*rho;
xstar=X-lag0(X,1)*rho;
ystar=ystar(2:end,:);
xstar=xstar(2:end,:);
M=inv(inv(Sigma0)+(1/sigma2)*(xstar'*xstar))*(inv(Sigma0)*B0+(1/sigma2)*xstar'*ystar);
V=inv(inv(Sigma0)+(1/sigma2)*(xstar'*xstar));
chck=-1;
while chck<0
B=M+(randn(1,3)*chol(V))';
b=[B(2) B(3);1   0];
ee=max(abs(eig(b)));
if ee<=1
    chck=1;
end
end

%step 3 compute rho
y=Y-X*B;
x=lag0(y,1);
y=y(2:end);
x=x(2:end);
MM=inv(inv(Sigma0r)+(1/sigma2)*(x'*x))*(inv(Sigma0r)*rho0+(1/sigma2)*x'*y);
VV=inv(inv(Sigma0r)+(1/sigma2)*(x'*x));
%draw rho but again ensure stationarity
chck=-1;
while chck<0
rho=MM+(randn(1,1)*chol(VV))';

ee=abs(rho);
if ee<=1
    chck=1;
end
end



%step 3 sample sigma2 conditional on B from IG(T1,D1);
%compute residuals
resids=ystar-xstar*B;
%compute posterior df and scale matrix
T1=T0+T;
D1=D0+resids'*resids;
%draw from IG
z0=randn(T1,1);
z0z0=z0'*z0;
sigma2=D1/z0z0;

if i>burn
    %compute forecast for 12 quarters
    yhat=zeros(14,1);
    vhat=zeros(14,1);
    yhat(1:2)=Y(end-1:end); %starting values
        cfactor=sqrt(sigma2);  %standard deviation of the shocks
    for m=3:14
        vhat(m)=vhat(m-1)*rho+randn(1,1)*cfactor;
        yhat(m)=[1 yhat(m-1) yhat(m-2)]*B+vhat(m);
        
    end
    %save
    out1=[out1 [Y;yhat(3:end)]];
    out2=[out2;B'];
    out3=[out3;rho];
    out4=[out4;sigma2];
end
end


figure(1)
TT=1947.25:0.25:2012.5;
out2x=prctile(out1',[10 20 30 40 50 60 70 80 90])';
plot(TT,out2x);

xlim([2000 2013])

save output.mat out2 out3 out4



