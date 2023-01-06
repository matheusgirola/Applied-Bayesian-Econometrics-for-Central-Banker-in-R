clear
addpath('functions');


% a bi-variate VAR with dummy variable implementation of priors
%load data
data=xlsread('\data\datain.xls'); %data for US GDP growth and inflation 1948q1 2010q4

N=size(data,2);
L=2;   %number of lags in the VAR
Y=data;
X=[ones(size(Y,1),1) lag0(data,1) lag0(data,2) ];
Y=Y(3:end,:);
X=X(3:end,:);
T=rows(X);

%compute standard deviation of each series residual via an ols regression
%to be used in setting the prior
%first variable
y=Y(:,1);
x=[ones(T,1) X(:,2)]; 
b0=inv(x'*x)*(x'*y);
s1=sqrt(((y-x*b0)'*(y-x*b0))/(rows(y)-2));  %std of residual standard error
%second variable
y=Y(:,2);
x=[ones(T,1) X(:,3)]; 
b0=inv(x'*x)*(x'*y);
s2=sqrt(((y-x*b0)'*(y-x*b0))/(rows(y)-2));  

%mean of the data
mu=mean(Y);



%specify parameters of the minnesota prior
tau=0.1;  %controls prior on own first lags
d=1;    %decay for higher lags
lamdac=1;  %prior for the constant
lamda=1;   %sum of coefficients unit roots
delta=1;  %cointegration prior



%specify dummy observations for first lag
yd1=[(1/tau)*s1 0;
    0        (1/tau)*s2];

xd1=[0 (1/tau)*s1 0 0 0;
    0 0 (1/tau)*s2 0 0];

%specify dummies for second lag

yd2=[0  0;
    0   0];

xd2=[0 0 0 (1/tau)*s1*2^d  0;
    0 0  0  0 (1/tau)*s2*2^d ];

%specify priors for the constants

yd3=[0  0;
    0   0];

xd3=[1/lamdac 0 0 0 0;
    1/lamdac 0  0  0 0];

%specify sum of coefficient dummies

yd4=[lamda*mu(1) 0;
     0   lamda*mu(2)];
 
 xd4=[0 lamda*mu(1) 0 lamda*mu(1) 0;
     0   0   lamda*mu(2) 0 lamda*mu(2)];
 
 %specify common stochastic trend dummies
 yd5=[delta*mu(1) delta*mu(2)];
 xd5=[delta delta*mu(1) delta*mu(2) delta*mu(1) delta*mu(2)]; 
 

 %specify dummy variables for covariance matrix
 yd6=[s1 0;
     0   s2];
 
 xd6=[0 0 0 0 0;
     0   0  0 0 0];
%all dummy observations

yd=[yd1;yd2;yd3;yd4;yd5;yd6];
xd=[xd1;xd2;xd3;xd4;xd5;xd6];

Ystar=[Y;yd];
Xstar=[X;xd];
Tstar=rows(Xstar);
%compute posterior mean
betahat=inv(Xstar'*Xstar)*(Xstar'*Ystar); 
%compute inital value of sigma
e=Ystar-Xstar*betahat;
sigma=(e'*e)/Tstar;


REPS=2000;
BURN=1000;



%gibbs algorithm
out1=[];
out2=[];
for i=1:REPS
    M=vec(betahat);
    V=kron(sigma,inv(Xstar'*Xstar));
    %draw beta
    beta=M+(randn(1,N*(N*L+1))*chol(V))';
    %draw sigma
    e=Ystar-Xstar*reshape(beta,N*L+1,N);
    scale=e'*e;
    sigma=IWPQ(Tstar,inv(scale));
    
    if i>BURN;
%forecast
yhat=zeros(14,2);
   yhat(1:2,:)=Y(end-1:end,:);
   for i=3:14
    yhat(i,:)=[1 yhat(i-1,:) yhat(i-2,:)]*reshape(beta,N*L+1,N)+randn(1,N)*chol(sigma);
end
out1=[out1 [Y(:,1);yhat(3:end,1)]];
out2=[out2 [Y(:,2);yhat(3:end,2)]];
end
end
    

TT=1948.75:0.25:2014;
subplot(1,2,1)
plot(TT,prctile(out1,[50 10 20 30 70 80 90],2))
xlim([1995 2015])
title('GDP Growth');


subplot(1,2,2)
plot(TT,prctile(out2,[50 10 20 30 70 80 90],2))
xlim([1995 2015])
legend('Median Forecast','10th percentile','20th percentile','30th percentile','70th percentile','80th percentile','90th percentile');
title('Inflation');




