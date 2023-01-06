clear;
addpath('functions');
%generate artificial data
T=500;
B1=0.2;
B2=0.9;
C1=1;
C2=-1;
S1=3;
S2=1;
P=[0.95 0.05;0.05 0.95];
strue=zeros(T,2);
strue(1,1)=1;
strue=simS(strue,P);
e=randn(T,1);
Y=zeros(T,1);
X=zeros(T,1);
for i=2:T;
    X(i,:)=Y(i-1,:);
    if strue(i,1)==1
    Y(i)=[X(i,:) 1]*[B1 C1]'+e(i)*sqrt(S1);
    else
    Y(i)=[X(i,:) 1]*[B2 C2]'+e(i)*sqrt(S2);
    end
end
%data
y=Y;
x=[X ones(T,1)];

%specify starting values
phi1=[0.5;1];   %regime 1 coefficients
phi2=[0.8;-1];   %regime 2 coefficients
sig1=3;       %regime 1 variance
sig2=1;       %regime 2 variance
p=0.95;
q=0.95;
pmat=[p 1-q;1-p q];
ncrit=10; %each regime should have ncrit obs
%set Priors

%coefficients
B0=zeros(2,1); %prior mean
Sigma0=eye(2); %prior variance

%variances
d0=0.1; %prior scale
v0=1;   %prior df

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
    iS1=1/sig1;
    iS2=1/sig2;
    lik=0;
    filter=zeros(T,2);
    for j=1:T
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
           out4=[out4;[p q]];
           count=count+1;
       end
    
   end
   igibbs=igibbs+1;
     disp(sprintf(' Replication %s , %s Saved Draws %s. ', ... 
             num2str(igibbs), num2str(count) ));

end
  
figure(1)
subplot(5,2,1);
hist(out1(:,1),50);
vline(B1)
title('Coefficient regime 1');
axis tight
subplot(5,2,2);
hist(out1(:,3),50);
title('Coefficient regime 2');
vline(B2)
axis tight
subplot(5,2,3);
hist(out1(:,2),50);
vline(C1)
title('Intercept regime 1');
axis tight
subplot(5,2,4);
hist(out1(:,4),50);
title('Intercept regime 2');
vline(C2)
axis tight
subplot(5,2,5);
hist(out2(:,1),50);
vline(S1)
title('\sigma_{1}');
axis tight
subplot(5,2,6);
hist(out2(:,2),50);
vline(S2)
title('\sigma_{2}');
axis tight
subplot(5,2,7);
hist(out4(:,1),50);
title('P_{00}');
vline(P(1,1))
axis tight
subplot(5,2,8);
hist(out4(:,2),50);
title('p_{11}');
vline(P(2,2))
axis tight
subplot(5,2,[9 10])
temp=mean(out3,1);
plot(temp,'c','LineWidth',2);
hold on
plot(strue(:,2),'k','LineWidth',2)
title('Probability of Regime 1');
legend('Estimate','True')
axis tight


