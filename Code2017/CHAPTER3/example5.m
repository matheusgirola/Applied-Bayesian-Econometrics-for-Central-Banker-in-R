clear;
addpath('functions')


%generate artificial data 
nobs=996; %996 months 332 quarters
btrue=[0.95 0.1;
       0.1 0.95;
       -0.1 0;
       0    -0.1;
       -0.05  0;
       0    -0.05;
       0     0];
   
   sigmatrue=[2  1;
              1 2];
      
 datatrue=zeros(nobs,2);
 for j=4:nobs
 datatrue(j,:)=[datatrue(j-1,:) datatrue(j-2,:) datatrue(j-3,:) 1]*btrue+randn(1,2)*chol(sigmatrue);
 end
 %assume first variable is subject to temporal aggregation
 dataQ=zeros(nobs/3,1); %quarterly data Y
 jj=1;
 for j=1:3:nobs
     tmp=datatrue(j:j+2,1);
     dataQ(jj,:)=mean(tmp);
     jj=jj+1;
 end
dataM=datatrue(:,2); %monthly data X

%arrange data
%put missing observations
dataN=[  nan(rows(dataQ),2) dataQ(:,1) ];  %puts NANs for missing obs
dataN=vecr(dataN);
data0=[ zeros(rows(dataQ),2) dataQ(:,1) ];  %same as above but zeros for missing
data0=vecr(data0);
%initial value of data just repeated observations
dataX=repmat(dataQ(:,1),1,3);
dataX=vecr(dataX); %
data=[dataX dataM];
dataid=[ dataN dataM];
dataid0=[ data0 dataM];
mid=isnan(dataid);  %id for missing obs


N=cols(data);
REPS=11000;
BURN=10500;

L=3;  %lags
Y=data;
X=prepare(data,L); %X=[Y(-1),Y(-2)...constant]
Y=Y(L+1:end,:);
X=X(L+1:end,:);
dataid0=dataid0(L+1:end,:);
dataM=dataM(L+1:end,:);



T=rows(X);

%initial values for VAR coefficients

b0=X\Y;  %ols
e0=Y-X*b0;

sigma=eye(N);

%priors for VAR coefficients (Banbura et.al)
lamdaP  = 1;
tauP    = 10*lamdaP;
epsilonP= 1;
muP=mean(Y)';
sigmaP=[];
deltaP=[];
e0=[];
for i=1:N
    ytemp=Y(:,i);
    xtemp=[lag0(ytemp,1) ones(rows(ytemp),1)];
    ytemp=ytemp(2:end,:);
    xtemp=xtemp(2:end,:);
    btemp=xtemp\ytemp;
    etemp=ytemp-xtemp*btemp;
    stemp=etemp'*etemp/rows(ytemp);
    if abs(btemp(1))>1
        btemp(1)=1;
    end
    deltaP=[deltaP;btemp(1)];
    sigmaP=[sigmaP;stemp];
    e0=[e0 etemp];
end

%dummy data to implement priors see http://ideas.repec.org/p/ecb/ecbwps/20080966.html
[yd,xd] = create_dummies(lamdaP,tauP,deltaP,epsilonP,L,muP,sigmaP,N);


%Initial values for the Kalman filter B0/0
beta0=[];
for j=0:L-1
    beta0=[beta0 Y(L-j,:)];
end
P00=eye(cols(beta0))*0.1;  %P[0/0]


% Gibbs sampler
gibbs1=1;
for gibbs=1:REPS



%step 1 Draw VAR coefficients  
X0=[X;xd]; %add dummy obs
Y0=[Y;yd];
mstar=vec(X0\Y0);
vstar=kron(sigma,invpd(X0'*X0));
chck=-1;
while chck<0
varcoef=mstar+(randn(1,N*(N*L+1))*chol(vstar))'; %draw but keep stable
ee=stability(varcoef,N,L);
if ee==0;
    chck=1;
end
end


%step 2 Draw VAR covariance
 resids=Y0-X0*reshape(varcoef,N*L+1,N);
scaleS=(resids'*resids);
sigma=iwpQ(T,invpd(scaleS)); %draw for inverse Wishart

%step 3 Carter Kohn algorithm  to draw monthly data

ns=cols(P00);
[F,MUx]=comp(varcoef,N,L,1); %companion form for coefficients
Q=zeros(ns,ns);
Q(1:N,1:N)=sigma; %companion form for covariance

%Carter and Kohn algorithm to draw the factor
beta_tt=zeros(T,ns);          %will hold the filtered state variable
ptt=zeros(T,ns,ns);    % will hold its variance
% %%%%%%%%%%%Step 6a run Kalman Filter
beta11=beta0;
p11=P00;
for i=1:T
nanid=mid(i,1); %checks if data on GDP is missing
if nanid==1 %missing
 H=[0 0 0 0 0 0;
    0 1 0  0  0  0];
  
    rr=zeros(1,N);
    rr(1)=1e10;  %big variance so missing data ignored
    R=diag(rr);
else  %valid  observation for first variable every 3rd month
     H=[1/3 0 1/3 0 1/3 0;
    0 1 0  0  0  0];
   
    rr=zeros(1,N);
    R=diag(rr);
  
  
end
    
x=H;




    %Prediction
beta10=MUx+beta11*F';
p10=F*p11*F'+Q;
yhat=(x*(beta10)')';                                               
eta=dataid0(i,:)-yhat;
feta=(x*p10*x')+R;
%updating
K=(p10*x')*invpd(feta);
beta11=(beta10'+K*eta')';
p11=p10-K*(x*p10);
ptt(i,:,:)=p11;
beta_tt(i,:)=beta11;
end
% Backward recursion to calculate the mean and variance of the distribution of the state
%vector
beta2 = zeros(T,ns);   %this will hold the draw of the state variable
bm2=beta2;
jv=1:2; %index of non singular block
jv1=[1 3 5]; %state variables to draw, 3, 5 are lagged states

wa=randn(T,ns);

i=T;  %period t
p00=squeeze(ptt(i,jv1,jv1)); 
beta2(i,:)=beta_tt(i,:);
beta2(i,jv1)=mvnrnd(beta_tt(i:i,jv1),p00,1);%beta_tt(i:i,jv1)+(wa(i:i,jv1)*cholx(p00));   %draw for beta in period t from N(beta_tt,ptt)
q=Q(jv,jv);
mu=MUx(jv);
f=F(jv,:);
%periods t-1..to .1
for i=T-1:-1:1
   
pt=squeeze(ptt(i,:,:));
bm=beta_tt(i:i,:)+(pt*f'*invpd(f*pt*f'+q)*(beta2(i+1:i+1,jv)-mu-beta_tt(i,:)*f')')';  
pm=pt-pt*f'*invpd(f*pt*f'+q)*f*pt;  
beta2(i,:)=bm;
beta2(i:i,jv1)=mvnrnd(bm(jv1),pm(jv1,jv1),1);     %bm(jv1)+(wa(i:i,jv1)*cholx(pm(jv1,jv1)));  
bm2(i,:)=bm;
end

out=beta2(:,1); %draw of monthly data








    
datax=[out dataM];
    Y=datax;
X=prepare(Y,L);
Y=Y(L+1:end,:);
X=X(L+1:end,:);



    
disp(sprintf('Iteration Number= %s ', num2str(gibbs)));

if gibbs>=BURN

 dmat(:,gibbs1)=out;
 bmat(gibbs1,:)=varcoef;
 smat(gibbs1,:,:)=sigma;
gibbs1=gibbs1+1;  
end

end

figure(1)
tmp=prctile(dmat,[50 ],2);
plot(tmp,'r');hold on;
plot(datatrue(4:end,1),'k')
legend('Posterior Median','True Data');

 