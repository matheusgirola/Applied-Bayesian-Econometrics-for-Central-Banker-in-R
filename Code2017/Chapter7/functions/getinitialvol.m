function [outh,outf,outg,outA]=getinitialvol(data0,REPS,BURN,T0,L)

N=cols(data0);
data=data0;

Y=data;
trend=(0:length(Y)-1)';

%take lags
X=[];
for j=1:L
X=[X lag0(data,j) ];
end
X=[X ones(rows(X),1) trend];
Y=Y(L+1:end,:);
X=X(L+1:end,:);
T=rows(X);

%initial training sample
Y0=Y(1:T0,:);
X0=X(1:T0,:);


%priors for VAR coefficients
lamdaP=0.1;  %This controls the tightness of the priors on the first lag
tauP=10*lamdaP;  % this controls the tightness of the priors on sum of coefficients
epsilonP=1/1000;  % this controls tightness of the prior on the constant
muP=mean(Y0)';
sigmaP=[];
deltaP=[];
e0=[];
for i=1:N
    ytemp=Y0(:,i);
    xtemp=[lag0(ytemp,1) ones(rows(ytemp),1)];
    ytemp=ytemp(2:end,:);
    xtemp=xtemp(2:end,:);
    btemp=xtemp\ytemp;
    etemp=ytemp-xtemp*btemp;
    stemp=etemp'*etemp/rows(ytemp);
    if btemp(1)>1
        btemp(1)=1;
    end
    deltaP=[deltaP;btemp(1)];
    sigmaP=[sigmaP;sqrt(stemp)];
    e0=[e0 etemp];
end
%dummy data to implement priors see http://ideas.repec.org/p/ecb/ecbwps/20080966.html
[yd,xd] = create_dummiestrendx(lamdaP,tauP,deltaP,epsilonP,L,muP,sigmaP,N);

%prior mean and variance using the dummy data
B0=xd\yd;
E0=yd-xd*B0;
 SIGMA0=kron(E0'*E0,pinv(xd'*xd)); %prior variance
% SIGMA0=kron(eye(N),pinv(xd'*xd)); %prior variance

B0=vec(B0); %prior mean


%priors for the elements of A matrix
s0=(e0'*e0)/T0;
C0=chol(s0);
C0=inv(C0./repmat(diag(C0),1,N))'; %prior mean
SC0=10;  %prior variance

%priors for the initial condition of the stochastic volatility
MU0=log(diag(s0)); %prior mean
SV0=1;           %prior variance

%Remove training sample. This is the data used for estimation
  Y=Y(T0+1:end,:);
  X=X(T0+1:end,:);
  T=rows(Y);
  %rough guess for stochastic volatility
 
hlast=(diff(Y).^2)+0.0001;
hlast=[hlast(1:2,:);hlast];  %rough intial guess for svol
errors=diff(Y); %initial guess for VAR residuals
errors=[errors(1,:);errors(1:end,:)];
g=ones(N,1);  %rough guess for the variance of the transition equation
g0=0.01^2;  %prior scale parameter for inverse gamma
Tg0=1;     %prior degrees of freedom
SV=zeros(REPS-BURN,T+1,N);
h01=hlast;
h0=hlast;
errors0=errors;
for i=1:50
    for j=1:N
    h01(:,j)=getsvol(h0(:,j),g(j),log(s0(j,j)),10,errors0(:,j));
    end
    h0=h01;
end
hlast=h0;
B0h=[0.7;0];
Sigma0h=diag([0.5 0.00005]);

MUSVOL=zeros(N,1);
FSVOL=zeros(N,1);

beta20=vec(X\Y);






tic;
jgibbs=1;
for igibbs=1:REPS
    igibbs
%step 1 of the Gibbs algorithm. Sample the A matrix  (See step 3 in the pdf
%file)

 amatx=[];

for j=2:N
    
%v2=-a1*v1+sqrt(exp(h2))*e2
ytemp=errors(:,j);  
xtemp=errors(:,1:j-1)*(-1); 
ytemp=ytemp./sqrt(hlast(2:rows(hlast),j)); %remove heteroscedasticity
xtemp=xtemp./repmat(sqrt(hlast(2:rows(hlast),j)),1,cols(xtemp));%remove heteroscedasticity
%prior means and variances
A0=C0(j,1:j-1)';
V00=diag(abs(A0))*SC0;

atemp=getregx(ytemp,xtemp,A0,V00,1); 

amatx=[amatx atemp'];
end
  
%step 2 of the algorithm: Sample stochastic volatilities  
A=chofac(N,amatx');
epsilon=errors*A';  %orthogonal residuals
%sample stochastic vol for each epsilon using the MH algorithm in Jaquier
%Polson and Rossi
  hnew=[];
  for i=1:N
      htemp=getsvolx(hlast(:,i),g(i),MU0(i),SV0,epsilon(:,i),MUSVOL(i),FSVOL(i)); 
      
      hnew=[hnew htemp];
  end
hlast=hnew;




%draw AR parameters
 gerrors=[];
 for jj=1:N
     ytemp=log(hlast(:,jj));
     xtemp=[lag0(ytemp,1) ones(rows(ytemp),1)];
     ytemp=ytemp(2:end,:);
     xtemp=xtemp(2:end,:);


MM=inv(inv(Sigma0h)+(1/g(jj))*(xtemp'*xtemp))*(inv(Sigma0h)*B0h+(1/g(jj))*xtemp'*ytemp); 
VV=inv(inv(Sigma0h)+(1/g(jj))*(xtemp'*xtemp));
chck=-1;
while chck<0                     %check for stability
BB=MM+(randn(1,2)*chol(VV))';
ee=max(abs(BB(1)));
if ee<=1
    chck=1;
end
end
     FSVOL(jj)=BB(1);
     MUSVOL(jj)=BB(2);
     gerrors=[gerrors ytemp-xtemp*BB];
 end













 
 %step 3 of the algorithm: Sample g from the IG distribution
 for i=1:N
    
g(i)=IG(Tg0,g0,gerrors(:,i));  %draw from the inverse Gamma distribution
 end

 %step 4: Sample VAR coefficients
%  [beta2,errors]=kcohnxx(Y,X,amatx,hlast,B0',SIGMA0,L);  %Conditional distribution of the VAR coefficients using the Kalman filter
 [beta2,errors,roots,problem]=carterkohnvar(Y,X,0,invpd(A),hlast,B0',SIGMA0,L,1,50,2);
problem

if problem
    beta2=beta20;
else
    beta20=beta2;
    
end
 if igibbs>BURN
     out1(jgibbs,:)=beta2;
    out2(jgibbs,1:T+1,1:N)=hlast(1:end,:);
    out4(jgibbs,1:T,1:(N*(N-1))/2)=repmat(amatx,T,1);
 
    out7(jgibbs,1:N)=g';
    out8(jgibbs,1:N*2)=[FSVOL' MUSVOL'];
     jgibbs=jgibbs+1;
 end
end

outh=squeeze(mean(out2,1));
outf=mean(out8);
outg=mean(out7);
tmp=squeeze(mean(out4,1));
outA=tmp(end,:);
