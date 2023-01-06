clear;
addpath('functions');
%generate artificial data from a MSTVTP model
T=500;
B1=0.2;
B2=0.9;
C1=1;
C2=-1;
S1=3;
S2=1;
GAMMA0=-1;
GAMMA1=1;
LAMBDA0=10;

strue=zeros(T,1);
strue(1,1)=1;

Z=getar(0.9,T);%randn(T,1);
e=randn(T,1);
Y=zeros(T,1);
X=zeros(T,1);
SSTAR=zeros(T,1);
ptrue=zeros(T,1);
qtrue=zeros(T,1);
for i=2:T;
    X(i,:)=Y(i-1,:);
    SSTAR(i,:)=GAMMA0+Z(i,:)*LAMBDA0+GAMMA1*strue(i-1,1)+randn(1,1);
    if SSTAR(i,:)>=0
        strue(i,1)=1;
   
    end
    %transition probabilities
    ptrue(i)=normcdf((-GAMMA0-Z(i,:)*LAMBDA0));
    qtrue(i)=1-normcdf((-GAMMA0-Z(i,:)*LAMBDA0-GAMMA1));
    
    
    if strue(i,1)==0
    Y(i)=[X(i,:) 1]*[B1 C1]'+e(i)*sqrt(S1);
    else
    Y(i)=[X(i,:) 1]*[B2 C2]'+e(i)*sqrt(S2);
    end
end
%data
y=Y;
x=[X ones(T,1)];
z=Z;

%specify starting values
phi1=[0.5;1];   %regime 1 coefficients
phi2=[0.8;-1];   %regime 2 coefficients
sig1=3;       %regime 1 variance
sig2=1;       %regime 2 variance
gamma=[-1 0 1]'; %coefficients of prob equation
pp=repmat(0.95,T,1);
qq=repmat(0.95,T,1);
ncrit=10; %each regime should have ncrit obs
%set Priors

%coefficients
B0=zeros(2,1); %prior mean
Sigma0=eye(2); %prior variance

%variances
d0=0.1; %prior scale
v0=1;   %prior df

%transition probabilities
GAMMA00=zeros(3,1); %prior mean coefficients of probability equation
SGAMMA0=eye(3).*1000;

out1=[];  %save coefficients
out2=[];  %save variances
out3=[];  %save S
out4=[]; %save p00
out5=[]; %save p11
out6=[]; %save gamma
out7=[]; %save sstar


REPS=20000;
BURN=15000;
igibbs=1;
count=1;
while count<REPS-BURN
    
 

%step 1: sample S[t]

 %%%%%%%%%%%%%%%%Run Hamilton Filter%%%%%%%%%%%%%%%%
   %unconditional probabilities
pmat=[pp(1) 1-qq(1);
      1-pp(1) qq(1)];
A = [(eye(2)-pmat);ones(1,2)];
           EN=[0;0;1];
           ett11= pinv(A'*A)*A'*EN;
    iS1=1/sig1;
    iS2=1/sig2;
    lik=0;
    filter=zeros(T,2);
    for j=1:T
        pmat=[pp(j) 1-qq(j); %TVP transition prob
      1-pp(j) qq(j)];
        em1=y(j)-x(j,:)*phi1; 
        em2=y(j)-x(j,:)*phi2; 
        neta1=(1/sqrt(sig1))*exp(-0.5*(em1*iS1*em1'));%F(Y\S=0)
        neta2=(1/sqrt(sig2))*exp(-0.5*(em2*iS2*em2'));%F(Y\S=1)
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
          pmat=[pp(t+1) 1-qq(t+1); %TVP transition prob
      1-pp(t+1) qq(t+1)];
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
 
 
 %step 2 sample the transition Probabilties
    
   %step 2a Sample sstar
sstar=zeros(T,1);
Slag=lag0(S,1);
Slag(1)=Slag(2);
zall=[ones(T,1) z Slag];
mm=zall*gamma;
for t=1:T
          
if S(t)==1
sstar(t)= normlt_rnd(mm(t),1,0);%draw from left truncated normal N(mm,1)
elseif S(t)==0;
sstar(t)= normrt_rnd(mm(t),1,0); %draw from right truncated normal N(mm,1)
end

end

  %step 2 b Calculate pp,qq
  pp=normcdf(-zall(:,1:end-1)*gamma(1:end-1));
  qq=1-normcdf(-zall(:,1:end-1)*gamma(1:end-1)-gamma(end));
  
  %step 2c Sample gamma
  yy=sstar;
  xx=zall;
  V=inv(inv(SGAMMA0)+(xx'*xx));
   M=V*(inv(SGAMMA0)*GAMMA00+xx'*yy); 
   gamma=M+(randn(1,3)*chol(V))';
    
    %step 3 sample beta
    % Select data in regime 1
    id=find(S==0);
    y1=y(id);
    x1=x(id,:);
    M=inv(inv(Sigma0)+(1/sig1)*(x1'*x1))*(inv(Sigma0)*B0+(1/sig1)*x1'*y1); 
    V=inv(inv(Sigma0)+(1/sig1)*(x1'*x1));
    phi1=M+(randn(1,2)*chol(V))';
    %Select data in regime 2
    id=find(S==1);
    y2=y(id);
    x2=x(id,:);
    M=inv(inv(Sigma0)+(1/sig2)*(x2'*x2))*(inv(Sigma0)*B0+(1/sig2)*x2'*y2); 
    V=inv(inv(Sigma0)+(1/sig2)*(x2'*x2));
    phi2=M+(randn(1,2)*chol(V))';
    
    
    %step 4 sample sigma
    
    %residuals regime 1
    e1=y1-x1*phi1;
    T1=v0+rows(e1);
    D1=d0+e1'*e1;
    %draw from IG
   z0=randn(T1,1);
    z0z0=z0'*z0;
   sig1=D1/z0z0;
   %residuals regime 2
    e2=y2-x2*phi2;
    T2=v0+rows(e2);
    D2=d0+e2'*e2;
    %draw from IG
   z0=randn(T2,1);
    z0z0=z0'*z0;
   sig2=D2/z0z0;
   
   
   
   %save and impose regime identification
   if igibbs>BURN
       chck=phi1(2,1)>phi2(2,1); %constant bigger in regime 1
       if chck
           out1=[out1;([phi1' phi2'])];
           out2=[out2;([sig1 sig2 ])];
           out3=[out3;S'];
           out4=[out4;pp'];
           out5=[out5;qq'];
           out6=[out6;gamma'];
           out7=[out7;sstar'];
           count=count+1;
       end
    
   end
   igibbs=igibbs+1;
     disp(sprintf(' Replication %s , %s Saved Draws %s. ', ... 
             num2str(igibbs), num2str(count) ));

end
  
figure(1)
subplot(8,2,1);
hist(out1(:,1),50);
vline(B1)
title('Coefficient regime 1');
axis tight
subplot(8,2,2);
hist(out1(:,3),50);
title('Coefficient regime 2');
vline(B2)
axis tight
subplot(8,2,3);
hist(out1(:,2),50);
vline(C1)
title('Intercept regime 1');
axis tight
subplot(8,2,4);
hist(out1(:,4),50);
title('Intercept regime 2');
vline(C2)
axis tight
subplot(8,2,5);
hist(out2(:,1),50);
vline(S1)
title('\sigma_{1}');
axis tight
subplot(8,2,6);
hist(out2(:,2),50);
vline(S2)
title('\sigma_{2}');
axis tight
subplot(8,2,7);
hist(out6(:,1),50);
title('\gamma_0');
vline(GAMMA0)
axis tight
subplot(8,2,8);
hist(out6(:,2),50);
title('\lambda');
vline(LAMBDA0)
axis tight
subplot(8,2,9);
hist(out6(:,3),50);
title('\gamma_1');
vline(GAMMA1)
axis tight

subplot(8,2,[11 12])
temp=mean(out3,1);
plot(temp,'c','LineWidth',2);
hold on
plot(strue(:,1),'k','LineWidth',2)
title('Probability of Regime 1');
legend('Estimate','True')
axis tight

subplot(8,2,[13 14])
temp=mean(out4,1);
plot(temp,'c','LineWidth',2);
hold on
plot(ptrue(:,1),'k','LineWidth',2)
title('p_00');
legend('Estimate','True')
axis tight

subplot(8,2,[15 16])
temp=mean(out5,1);
plot(temp,'c','LineWidth',2);
hold on
plot(qtrue(:,1),'k','LineWidth',2)
title('p_11');
legend('Estimate','True')
axis tight


