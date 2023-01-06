clear
addpath('functions');
data=xlsread('\data\datain.xls'); %data for US GDP growth and inflation 1948q1 2010q4
N=cols(data);
horizon=3;
path=[1;1;1]; %constrained values for X

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


B=X\Y;  %ols estimate
res=Y-X*B;
sigma=(res'*res)/T;
A0=chol(sigma);
%calculate impulse responses to be used to construct R
S=zeros(1,N);
S(1)=1; %shock to first eq
Z1=irfsim(B,N,L,A0,S,horizon+L);
S=zeros(1,N);
S(2)=1; %shock to 2nd eq
Z2=irfsim(B,N,L,A0,S,horizon+L);
%calculate unconditional forecast to be used to construct r
yhat1=zeros(horizon+L,N);
yhat1(1:L,:)=Y(end-L+1:end,:);
for i=L+1:horizon+L   
    x=[];
    for j=1:L
    x=[x yhat1(i-j,:)];
    end
    yhat1(i,:)=[x 1]*B;  
end
yhat1=yhat1(L+1:end,:);
%construct the R matrix
R=[Z1(1,2) Z2(1,2) 0  0  0  0;
    Z1(2,2) Z2(2,2) Z1(1,2) Z2(1,2) 0  0;
    Z1(3,2) Z2(3,2) Z1(2,2) Z2(2,2) Z1(1,2) Z2(1,2)];
%construct the r matrix
r=path-yhat1(:,2);
%compute the restricted structural shocks
ehat=R'*pinv(R*R')*r;
ehat=reshape(ehat,N,horizon)';
%compute the conditional forecast
yhat2=zeros(horizon+L,N);
yhat2(1:L,:)=Y(end-L+1:end,:);
for i=L+1:horizon+L   
    x=[];
    for j=1:L
    x=[x yhat2(i-j,:)];
    end
    yhat2(i,:)=[x 1]*B+ehat(i-L,:)*A0;  
end
yhat2=yhat2(L+1:end,:);

%Gibbs sampling algorithm
REPS=5000;
BURN=3000;
out1=[]; %will hold forecast for GDP
out2=[]; %will hold forecast for inflation
yhatg=yhat2; %initialise conditional forecast
sig=sigma; %initialise error covariance
for igibbs=1:REPS
    
    %step 1 DRAW VAR parameters
 datag=[data;yhatg]; %appended data
 YSTAR=datag;
    %take lags
XSTAR=[];
for j=1:L
XSTAR=[XSTAR lag0(datag,j) ];
end
XSTAR=[XSTAR ones(rows(XSTAR),1)];
YSTAR=YSTAR(L+1:end,:);
XSTAR=XSTAR(L+1:end,:);
T=rows(XSTAR);
%conditional mean
M=vec(XSTAR\YSTAR);
%conditional variance
V=kron(sig,inv(XSTAR'*XSTAR));
bg=M+(randn(1,N*(N*L+1))*chol(V))';
bg1=reshape(bg,N*L+1,N);
%draw sigma from the IW distribution
e=YSTAR-XSTAR*bg1;
scale=e'*e;
sig=IWPQ(T,inv(scale));
%A0 matrix
A0g=chol(sig);


%step 2 Construct conditional forecast 
%impulse responses
S=zeros(1,N);
S(1)=1; %shock to first eq
Z1=irfsim(bg1,N,L,A0g,S,horizon+L);
S=zeros(1,N);
S(2)=1; %shock to 2nd eq
Z2=irfsim(bg1,N,L,A0g,S,horizon+L);
%calculate unconditional forecast to be used to construct r
yhat1=zeros(horizon+L,N);
yhat1(1:L,:)=Y(end-L+1:end,:);
for i=L+1:horizon+L   
    x=[];
    for j=1:L
    x=[x yhat1(i-j,:)];
    end
    yhat1(i,:)=[x 1]*bg1;  
end
yhat1=yhat1(L+1:end,:);
%construct the R matrix
R=[Z1(1,2) Z2(1,2) 0  0  0  0;
    Z1(2,2) Z2(2,2) Z1(1,2) Z2(1,2) 0  0;
    Z1(3,2) Z2(3,2) Z1(2,2) Z2(2,2) Z1(1,2) Z2(1,2)];
%construct the r matrix
r=path-yhat1(:,2);

%compute the mean of the distribution of restricted structural shocks
MBAR=R'*pinv(R*R')*r;
%compute the variance of the distribution of restricted structural shocks
VBAR=R'*pinv(R*R')*R;
VBAR=eye(cols(VBAR))-VBAR;
%draw structural shocks from the N(MBAR,VBAR) distribution
 edraw=MBAR+(randn(1,rows(MBAR))*real(sqrtm(VBAR)))';
 edraw=reshape(edraw,N,horizon)';
%conditional forecast using new draw of shocks
yhatg=zeros(horizon+L,N);
yhatg(1:L,:)=Y(end-L+1:end,:);
for i=L+1:horizon+L   
    x=[];
    for j=1:L
    x=[x yhatg(i-j,:)];
    end
    yhatg(i,:)=[x 1]*bg1+edraw(i-L,:)*A0g;  
end
yhatg=yhatg(L+1:end,:);

if igibbs>BURN
    out1=[out1;[Y(:,1);yhatg(:,1)]'];
    out2=[out2;[Y(:,2);yhatg(:,2)]'];
end
end
  
TT=1948.75:0.25:2011+(.75);
subplot(1,2,1)
plot(TT,prctile(out1,[50 10 20 30 70 80 90],1))
xlim([1995 max(TT)+0.25])
title('GDP Growth');


subplot(1,2,2)
plot(TT,prctile(out2,[50 10 20 30 70 80 90],1))
xlim([1995 max(TT)+0.25])
legend('Median Forecast','10th percentile','20th percentile','30th percentile','70th percentile','80th percentile','90th percentile');
title('Inflation');
