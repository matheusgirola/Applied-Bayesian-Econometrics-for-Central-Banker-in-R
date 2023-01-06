clear
addpath('functions');
% a bi-variate VAR with a Minnesota Prior and Gibbs Sampling
%load data
data=xlsread('\data\datain.xls'); %data for US GDP growth and inflation 1948q1 2010q4
N=size(data,2);
L=2;   %number of lags in the VAR
Y=data;
X=[lag0(data,1) lag0(data,2) ones(rows(data),1) ]; 
Y=Y(3:end,:);
X=X(3:end,:);
T=rows(X);

%compute standard deviation of each series residual via an ols regression
%to be used in setting the prior
%first variable
y=Y(:,1);
x=[ones(T,1) X(:,1)]; 
b0=inv(x'*x)*(x'*y);
s1=sqrt(((y-x*b0)'*(y-x*b0))/(rows(y)-2));  %std of residual standard error
%second variable
y=Y(:,2);
x=[ones(T,1) X(:,2)]; 
b0=inv(x'*x)*(x'*y);
s2=sqrt(((y-x*b0)'*(y-x*b0))/(rows(y)-2));  


%specify parameters of the minnesota prior
lamda1=1; %controls the prior on own lags
lamda2=1;
lamda3=1;
lamda4=1;

%specify the prior mean of the coefficients of the Two equations of the VAR
B01=[1;0;0;0];
B02=[0;1;0;0];

B0=[B01;B02];

%Specify the prior variance of vec(B)
H=zeros(8,8);
%for equation 1  of the VAR
H(1,1)=(lamda1)^2;     %own lag
H(2,2)=((s1*lamda1*lamda2)/s2)^2;  %lag of other variable
H(3,3)=(lamda1/(2^lamda3))^2;    %own second lag
H(4,4)=((s1*lamda1*lamda2)/(s2*(2^lamda3)))^2;  %lag of other variable
%for equation 2 of the VAR
H(5,5)=((s2*lamda1*lamda2)/s1)^2;  %lag of other variable
H(6,6)=(lamda1)^2;     %own lag
H(7,7)=((s2*lamda1*lamda2)/(s1*(2^lamda3)))^2;  %lag of other variable
H(8,8)=(lamda1/(2^lamda3))^2;    %own second lag

%prior scale matrix for sigma the VAR covariance
S=eye(N);
%prior degrees of freedom
alpha=N+1;

%set priors for the long run mean which is a N by 1 vector
M0=[1 1]; %prior mean  
V0=eye(N)*0.001;  %prior variance

%starting values via OLS
betaols=inv(X'*X)*(X'*Y);
F=[betaols(1:N*L,:)';eye(N*(L-1),N*L)]; %companion form
C=zeros(rows(F),1);
C(1:N)=betaols(N*L+1,:)';
MU=inv(eye(rows(F))-F)*C;  %ols estimate of the mean inv(I-B)C
e=Y-X*betaols;
Sigma=(e'*e)/T;

Reps=10000;
burn=5000;
out1=[]; %will store forecast of GDP growth
out2=[]; %will store forecast of inflation
i=1;
for j=1:Reps
%demean the data
  Y0=data-repmat(MU(1:N)',rows(data),1);
  X0=[];
  for jj=1:L
  X0=[X0 lag0(Y0,jj) ];
  end
  Y0=Y0(L+1:end,:);
  X0=X0(L+1:end,:);
%step 1 draw the VAR coefficients
bols=vec(inv(X0'*X0)*(X0'*Y0));
M=inv(inv(H)+kron(inv(Sigma),X0'*X0))*(inv(H)*B0+kron(inv(Sigma),X0'*X0)*bols);
V=inv(inv(H)+kron(inv(Sigma),X0'*X0));
beta=M+(randn(1,N*(N*L))*chol(V))';
beta1=reshape(beta,N*L,N);
%draw sigma from the IW distribution
e=Y0-X0*beta1;
%scale matrix
scale=e'*e+S;
Sigma=IWPQ(T+alpha,inv(scale));

%step 3 draw MU the long run mean conditional on beta and sigma (see
  %Appendix A in villani.
  
  Y1=Y-X(:,1:end-1)*beta1;  
  U=eye(N);
  jj=1;
  for jx=1:L
      betai=beta1(jj:jj+N-1,:);
      U=[U;betai'];
      jj=jj+N;
  end
  D=[ones(T,1) -ones(T,L)];
  vstar1=inv(U'*kron(D'*D,inv(Sigma))*U+inv(V0)); %posterior variance
  mstar1=vstar1*(U'*vec(inv(Sigma)*Y1'*D)+inv(V0)*M0'); %posterior mean
  MU=mstar1+(randn(1,N)*chol(vstar1))';  %draw MU
  


if j>burn
    
    %forecast GDP growth and inflation for 3 years
  F=[beta1(1:N*L,:)';eye(N*(L-1),N*L)]; %companion form
  mu=[];
  for i=1:L
  mu=[mu;MU];
  end
  C=(eye(rows(F))-F)*mu;  %implied constant
    
   
   yhat=zeros(44,2);
   yhat(1:2,:)=Y(end-1:end,:);
   for i=3:44
    yhat(i,:)=C(1:N)'+[yhat(i-1,:) yhat(i-2,:)]*reshape(beta,N*L,N)+randn(1,N)*chol(Sigma);
end
out1=[out1 [Y(:,1);yhat(3:end,1)]];
out2=[out2 [Y(:,2);yhat(3:end,2)]];
end

end

TT=1948.75:0.25:2021.5;
subplot(1,2,1)
plot(TT,[mean(out1,2) prctile(out1,[50 10 20 30 70 80 90],2)])
xlim([1995 2022])
title('GDP Growth');


subplot(1,2,2)
plot(TT,[ mean(out2,2) prctile(out2,[50 10 20 30 70 80 90],2)])
xlim([1995 2022])
legend('Mean Forecast','Median Forecast','10th percentile','20th percentile','30th percentile','70th percentile','80th percentile','90th percentile');
title('Inflation');




