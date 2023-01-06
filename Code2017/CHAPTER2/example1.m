clear
addpath('functions');
% a bi-variate VAR with a Minnesota Prior and Gibbs Sampling
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
x=X(:,1:2); 
b0=inv(x'*x)*(x'*y);
s1=sqrt(((y-x*b0)'*(y-x*b0))/(rows(y)-2));  %std of residual standard error
%second variable
y=Y(:,2);
x=X(:,[1 3]); 
b0=inv(x'*x)*(x'*y);
s2=sqrt(((y-x*b0)'*(y-x*b0))/(rows(y)-2));  


%specify parameters of the minnesota prior
lamda1=1; %controls the prior on own lags
lamda2=1;
lamda3=1;
lamda4=1;

%specify the prior mean of the coefficients of the Two equations of the VAR
B01=[0;1;0;0;0];
B02=[0;0;1;0;0];

B0=[B01;B02];

%Specify the prior variance of vec(B)
H=zeros(10,10);
%for equation 1  of the VAR
H(1,1)=(s1*lamda4)^2;  %constant 
H(2,2)=(lamda1)^2;     %own lag
H(3,3)=((s1*lamda1*lamda2)/s2)^2;  %lag of other variable
H(4,4)=(lamda1/(2^lamda3))^2;    %own second lag
H(5,5)=((s1*lamda1*lamda2)/(s2*(2^lamda3)))^2;  %lag of other variable
%for equation 2 of the VAR
H(6,6)=(s2*lamda4)^2;  %constant
H(7,7)=((s2*lamda1*lamda2)/s1)^2;  %lag of other variable
H(8,8)=(lamda1)^2;     %own lag
H(9,9)=((s2*lamda1*lamda2)/(s1*(2^lamda3)))^2;  %lag of other variable
H(10,10)=(lamda1/(2^lamda3))^2;    %own second lag

%prior scale matrix for sigma the VAR covariance
S=eye(N);
%prior degrees of freedom
alpha=N+1;

%starting values for the Gibbs sampling algorithm
Sigma=eye(N);
betaols=vec(inv(X'*X)*(X'*Y));

Reps=10000;
burn=5000;
out1=[]; %will store forecast of GDP growth
out2=[]; %will store forecast of inflation
i=1;
for j=1:Reps

%step 1 draw the VAR coefficients
M=inv(inv(H)+kron(inv(Sigma),X'*X))*(inv(H)*B0+kron(inv(Sigma),X'*X)*betaols);
V=inv(inv(H)+kron(inv(Sigma),X'*X));
beta=M+(randn(1,N*(N*L+1))*chol(V))';

%draw sigma from the IW distribution
e=Y-X*reshape(beta,N*L+1,N);
%scale matrix
scale=e'*e+S;
Sigma=IWPQ(T+alpha,inv(scale));

if j>burn
    %forecast GDP growth and inflation for 3 years
   yhat=zeros(14,2);
   yhat(1:2,:)=Y(end-1:end,:);
   for i=3:14
    yhat(i,:)=[1 yhat(i-1,:) yhat(i-2,:)]*reshape(beta,N*L+1,N)+randn(1,N)*chol(Sigma);
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




