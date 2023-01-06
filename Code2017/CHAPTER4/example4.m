clear;
addpath('functions');
%generate artificial data
T=500;
N=2; %2 variables in VAR
L=1; % 1 Lag
B1=[0.2 -0.1 -1; 0.5 -0.1 -1];
B2=[0.5 0.1 1; 0.7 0.1 1];
S1=[3 -0.5;-0.5 3];
S2=[1 0.1;0.1 1];
P=[0.95 0.05;0.05 0.95];
strue=zeros(T,2);
strue(1,1)=1;
strue=simS(strue,P);
e=randn(T,N);
Y=zeros(T,N);
X=zeros(T,N*L+1);
for i=2:T;
    X(i,:)=[Y(i-1,:) 1];
    if strue(i,1)==1
    Y(i,:)=X(i,:)*B1'+e(i,:)*chol(S1);
    else
    Y(i,:)=X(i,:)*B2'+e(i,:)*chol(S2);
    end
end
%data
y=Y;
x=X;

%specify starting values
maxtrys=1000; %number of trys for stable draw
phiols=x\y;
phi1=vec(phiols);   %regime 1 coefficients
phi2=vec(phiols);   %regime 2 coefficients
phi10=phi1;
phi20=phi2;
sig1=eye(N)*3;       %regime 1 variance
sig2=eye(N);       %regime 2 variance
p=0.95;
q=0.95;
pmat=[p 1-q;1-p q];
ncrit=10; %each regime should have ncrit obs
%set Priors

% VAR coefficients and variance priors via dummy observations
lamdaP  = 10;
tauP    = 10*lamdaP;
epsilonP= 1/10000;
muP=mean(y)';
sigmaP=[];
deltaP=[];
e0=[];
for i=1:N
    ytemp=y(:,i);
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

%transition probabilities
u00=25; %p00~D(u11,u22)
u01=5;
u11=25; %p11~D(u22,u21)
u10=5;

out1=[];  %save coefficients
out2=[];  %save variances
out3=[];  %save S
out4=[]; %save p
REPS=10000;
BURN=5000;
igibbs=1;
count=1;
while count<REPS-BURN

%step 1: sample S[t]

 %%%%%%%%%%%%%%%%Run Hamilton Filter%%%%%%%%%%%%%%%%
   %unconditional probabilities

A = [(eye(2)-pmat);ones(1,2)];
           EN=[0;0;1];
           ett11= pinv(A'*A)*A'*EN;
    iS1=inv(sig1);
    iS2=inv(sig2);
    lik=0;
    filter=zeros(T,2);
    for j=1:T
        em1=y(j,:)-x(j,:)*reshape(phi1,N*L+1,N); 
        em2=y(j,:)-x(j,:)*reshape(phi2,N*L+1,N); 
        neta1=(1/sqrt(det(sig1)))*exp(-0.5*(em1*iS1*em1'));%F(Y\S=0)
        neta2=(1/sqrt(det(sig2)))*exp(-0.5*(em2*iS2*em2'));%F(Y\S=1)
        %%%Prediction Step%%%%
        ett10=pmat*ett11;
        %%%%Update Step%%%%
        ett11=ett10.*[neta1;neta2]; %joint density F(Y,S)
        fit=sum(ett11);           %Marginal density F(Y)
        ett11=(ett11)/fit;    %conditional density F(S\Y) the weights of the likelihood
        filter(j,1:2)=ett11';      %save filter probability ett  
        lik=lik+log(fit);      %save log likelihood
        
    end  
   
   
   
 check=-1;
 while check<0
   %backward recursion to sample from H(S[t]\S[t+1],y)
   S=zeros(T,1);
   %time T
   p1=filter(T,1);
   p2=filter(T,2);
   p=p1/(p1+p2);
   u=rand(1,1);
   S(T,1)=(u>=p);
  
   for t=T-1:-1:1
   if S(t+1)==0
p00=pmat(1,1)*filter(t,1);
p01=pmat(1,2)*filter(t,2);
elseif S(t+1)==1
p00=pmat(2,1)*filter(t,1);
p01=pmat(2,2)*filter(t,2);
   end
  u=rand(1,1);
  p=p00/(p00+p01);
  if u<p
      S(t)=0;
  else
      S(t)=1;
  end
   end
   
if sum(S==0)>=ncrit && sum(S==1)>=ncrit
    check=1;
end
 end
 
 
 %step 2 sample the transition matrix P
    
    tranmat=switchg(S+1,[1;2]); %calculate the number of regime switches
    N00=tranmat(1,1); %S(t-1)=0 S(t)=0
    N01=tranmat(1,2); %S(t-1)=0 S(t)=1
    N10=tranmat(2,1); %S(t-1)=1 S(t)=0
    N11=tranmat(2,2); %S(t-1)=1 S(t)=1
    %draw from the dirichlet density
    p0=drchrnd([N00+u00;N01+u01]);
    p=p0(1,1); %p00
    p0=drchrnd([N10+u10;N11+u11]);
    q=p0(2,1); %p11
    pmat=[p 1-q;1-p q]; %transition prob matrix

    
    %step 3 sample parameters
    % VAR parameters in regime 0
    id=find(S==0);
    Y1=y(id,:);
    X1=x(id,:);
    Y0=[Y1;yd];
    X0=[X1;xd];
  %conditional mean of the VAR coefficients
  mstar1=vec(X0\Y0);  %ols on the appended data
  xx=X0'*X0;
  ixx1=xx\eye(cols(xx));
   [ phi1,PROBLEM1] = getcoef( mstar1,sig1,ixx1,maxtrys,N,L ); %draw VAR coefficients
     if PROBLEM1
         phi1=phi01;
     else
         phi01=phi1;
     end
     
     %draw covariance
     e=Y0-X0*reshape(phi1,N*L+1,N);
    scale=e'*e;
    sig1=iwpQ(rows(Y0),inv(scale));  
    
    % VAR parameters in regime 1
    id=find(S==1);
    Y2=y(id,:);
    X2=x(id,:);
    Y0=[Y2;yd];
    X0=[X2;xd];
  %conditional mean of the VAR coefficients
  mstar2=vec(X0\Y0);  %ols on the appended data
  xx=X0'*X0;
  ixx2=xx\eye(cols(xx));
   [ phi2,PROBLEM2] = getcoef( mstar2,sig2,ixx2,maxtrys,N,L );
     if PROBLEM2
         phi2=phi02;
     else
         phi02=phi2;
     end
     
     %draw covariance
     e=Y0-X0*reshape(phi2,N*L+1,N);
    scale=e'*e;
    sig2=iwpQ(rows(Y0),inv(scale));  
   
   
   %save and impose regime identification
   if igibbs>BURN
       chck=log(det(sig1))>log(det(sig2)); %Total bigger in regime 0
       if chck
           out1(count,:)=[phi1' phi2'];
           out2(count,:,:)=[sig1 sig2 ];
           out3=[out3;S'];
           out4=[out4;[p q]];
           count=count+1;
       end
    
   end
   igibbs=igibbs+1;
     disp(sprintf(' Replication %s , %s Saved Draws %s. ', ... 
             num2str(igibbs), num2str(count) ));

end
  
figure(1)

temp=mean(out3,1);
plot(temp,'c','LineWidth',2);
hold on
plot(strue(:,2),'k','LineWidth',2)
title('Probability of Regime 1');
legend('Estimate','True')
axis tight


figure(2)
tmp=[vec(B1');vec(B2')];
for j=1:rows(tmp);
    subplot(6,2,j);
    hist(out1(:,j))
    vline(tmp(j));
    title(strcat('Coefficient:',num2str(j)));
end