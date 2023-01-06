clear
addpath('functions');
[ data0 junk ]=xlsread('.\data\datain.xls');
[ junk names ]=xlsread('.\data\names.xls');
%names=names(1,2:end);
index=xlsread('.\data\index.xls');
dindex=index(:,1); %dindex=1 for series that are log differenced dindex=3 differencing without logs
index=index(:,2);  %index=1 for 'fast moving' series


%first difference the data where appropriate
data=[];
for i=1:cols(data0);
    if dindex(i)==1
        dat=log(data0(:,i));
        dat=diff(dat)*100;
    elseif dindex(i)==3
        dat=diff(data0(:,i));
    else
        dat=data0(2:end,i);
    end
    data=[data dat];
end
%standardise the data
data=standardise(data);


%load policy rate and standardize it
z=xlsread('.\data\baserate.xls');
z=z(2:end);
z=standardise(z);



KK=3;  %number of factors
L=2;  %number of lags in the VAR
N=KK+1; %number of Variables in var K factors plus the interest rate
NN=cols(data);% size of the panel
T=rows(data)
%step 1 of the algorithm set starting values and priors

%get an intial guess for the factor via principal components
pmat=extract(data,KK);
beta0=[pmat(1,:) z(1) zeros(1,N)];  %state vector S[t-1/t-1]
ns=cols(beta0);
P00=eye(ns);  %P[t-1/t-1]
rmat=ones(NN,1); %arbitrary starting value for the variance of the idiosyncratic component
Sigma=eye(N);  %arbitrary starting value for the variance of VAR errors

%flat prior for the factor loadings,variances and VAR
reps=5000;
burn=4000;
mm=1;
for m=1:reps;

%gibbs sampling 


%step 2 sample factor loadings
fload=[];
floadr=[];
error=[];
for i=1:NN
    y=data(:,i);
    if index(i)==0
        x=pmat;
    else
        x=[pmat z];
    end
    M=inv(x'*x)*(x'*y);
    V=rmat(i)*inv(x'*x);
    %draw
    ff=M+(randn(1,cols(x))*cholx(V))';
    
    %save
    if index(i)==0;
        fload=[fload;ff'];
        floadr=[floadr;0];
    else
           fload=[fload;ff(1:end-1)'];
        floadr=[floadr;ff(end)];
    end
    error=[error y-x*ff];
end


%for identification top K by K block of fload is identity
fload(1:KK,1:KK)=eye(KK);
%for identification top K by 1 block of Floadr is zero
floadr(24:24+KK-1,1)=zeros(KK,1);

%step 3 sample variance of the idiosyncratic components from inverse
%wishart

rmat=[];
for i=1:NN
    rmati= IG(0,0,error(:,i));
    rmat=[rmat rmati];
end

%step 4 sample VAR coefficients
Y=[pmat z];
X=[lag0(Y,1) lag0(Y,2) ones(rows(Y),1)];
Y=Y(2:end,:);
X=X(2:end,:);

M=vec(inv(X'*X)*(X'*Y));  %conditional mean
V=kron(Sigma,inv(X'*X)); %conditional variance
chck=-1;                 %make sure VAR is stationary
while chck<0
beta=M+(randn(1,N*(N*L+1))*cholx(V))';  %draw for VAR coefficients
S=stability(beta,N,L);
if S==0
    chck=10;
end
end
beta1=reshape(beta,N*L+1,N);

errorsv=Y-X*beta1;

%sample VAR covariance
scale=errorsv'*errorsv;
Sigma=iwpQ(T,inv(scale));



%step 5 prepare matrices for the state space
%Y=H*factors+e
%factors=MU+F*Factors(-1)+v
%e~N(0,R)
%v~N(0,Q)

%matrix of factor loadings
H=zeros(NN,(KK+1)*L);
H(1:rows(fload),1:KK+1)=[fload floadr];
H(rows(floadr)+1,KK+1)=1;
%matrix R
R=diag([rmat 0]);
%vector MU
MU=[beta1(end,:)';zeros(N*(L-1),1)]';
%matrix F
F=[beta1(1:N*L,:)';eye(N*(L-1),N*L)];
%matrix Q
Q=zeros(rows(F),rows(F));
Q(1:N,1:N)=Sigma;



%Carter and Kohn algorithm to draw the factor
beta_tt=[];          %will hold the filtered state variable
ptt=zeros(T,ns,ns);    % will hold its variance
% %%%%%%%%%%%Step 6a run Kalman Filter
i=1;
x=H;
%Prediction
beta10=MU+beta0*F';
p10=F*P00*F'+Q;
yhat=(x*(beta10)')';                                                
eta=[data(i,:) z(i,:)]-yhat;
feta=(x*p10*x')+R;
%updating
K=(p10*x')*inv(feta);
beta11=(beta10'+K*eta')';
p11=p10-K*(x*p10);
beta_tt=[beta_tt;beta11];
ptt(i,:,:)=p11;
for i=2:T
    %Prediction
beta10=MU+beta11*F';
p10=F*p11*F'+Q;
yhat=(x*(beta10)')';                                               
eta=[data(i,:) z(i,:)]-yhat;
feta=(x*p10*x')+R;
%updating
K=(p10*x')*inv(feta);
beta11=(beta10'+K*eta')';
p11=p10-K*(x*p10);
ptt(i,:,:)=p11;
beta_tt=[beta_tt;beta11];
end
% Backward recursion to calculate the mean and variance of the distribution of the state
%vector
beta2 = zeros(T,ns);   %this will hold the draw of the state variable
jv1=1:KK; %index of state variables to extract
jv=1:N;
wa=randn(T,ns);
f=F(jv,:);
q=Q(jv,jv);
mu=MU(jv);
i=T;  %period t
p00=squeeze(ptt(i,jv1,jv1)); 
beta2(i,:)=beta_tt(i,:);
beta2(i,jv1)=beta_tt(i:i,jv1)+(wa(i:i,jv1)*cholx(p00));   %draw for beta in period t from N(beta_tt,ptt)
%periods t-1..to .1
for i=T-1:-1:1
pt=squeeze(ptt(i,:,:));
bm=beta_tt(i:i,:)+(pt*f'*inv(f*pt*f'+q)*(beta2(i+1:i+1,jv)-mu-beta_tt(i,:)*f')')';  
pm=pt-pt*f'*inv(f*pt*f'+q)*f*pt;  
beta2(i,:)=bm;

beta2(i:i,jv1)=bm(jv1)+(wa(i:i,jv1)*cholx(pm(jv1,jv1)));  

end


pmat=beta2(:,1:3);   %update the factors

if m>burn
    %compute impulse response
    A0=cholx(Sigma);
    yhat=zeros(36,N);
vhat=zeros(36,N);
vhat(3,1:N)=[0 0 0 1];

for i=3:36
 yhat(i,:)=[yhat(i-1,:) yhat(i-2,:) 1]*[beta1(1:N*L,:);zeros(1,N)]+vhat(i,:)*A0;
end

yhat1=yhat*H(:,1:KK+1)';  %impulse response for the panel

irfmat(mm,1:36,1:NN+1)=(yhat1);

mm=mm+1;
end
    
m
end



irf=prctile(irfmat,[50 16 84],1);

names{41}='Interest Rate';
figure(1)
j=1
for i=1:size(irf,3)
subplot(5,10,j)
plotx1(squeeze(irf(:,3:end,i))');
title(strcat('\fontsize{8}', names(i)))
j=j+1
axis tight
end