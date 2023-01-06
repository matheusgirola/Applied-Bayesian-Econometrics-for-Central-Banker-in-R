clear;
addpath('functions');
%generate artificial data
T=500;
B1=0.25;
B2=0.9;
MU1=3;
MU2=-1;
S1=1;
S2=2;
P0=[0.98 0.02;0.02 0.98];
P=matf(P0,2,1);
strue=zeros(T,4);
strue(1,1)=1;
strue=simS(strue,P);
e=randn(T,1);
Y=zeros(T,1);
X=zeros(T,1);
for i=2:T;
    X(i,:)=Y(i-1,:);
    if strue(i,1)==1
    Y(i)=(X(i,:)-MU1)*B1+MU1+e(i)*sqrt(S1);
    elseif strue(i,2)==1
    Y(i)=(X(i,:)-MU1)*B1+MU2+e(i)*sqrt(S2); 
    elseif strue(i,3)==1
    Y(i)=(X(i,:)-MU2)*B1+MU1+e(i)*sqrt(S1); 
    else
     Y(i)=(X(i,:)-MU2)*B1+MU2+e(i)*sqrt(S2);         
   
    end
end

y=Y;
x=X;

%specify starting values
phi=0.5; %AR coefficient
mu1=1;   %mean regime 0
mu2=0.1; %mean regime 1
sig1=1;       %regime 0 variance
sig2=1;       %regime 1 variance
p=0.95;
q=0.95;
pmat0=[p 1-q;1-p q];
pmat=matf(pmat0,2,1); %matf(P*,number of states, number of lagged states)
ncrit=10; %min number of obs in each regime

%set Priors

%AR coefficients
B0=zeros(1,1); %prior mean
Sigma0=eye(1); %prior variance

%mean
M0=zeros(2,1);
Sigma0M=eye(2)*10;

%variances
d0=0.1;
v0=1;

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

   % run hamilton filter 
 %unconditional probabilities
A = [(eye(4)-pmat);ones(1,4)];
EN=[0;0;0;0;1];
ett11= pinv(A'*A)*A'*EN;
isig1=1/sig1;
isig2=1/sig2;
isig3=1/sig1;
isig4=1/sig2;
lik=0;
filter=zeros(T,4);
for j=1:T
em1=(y(j)-mu1)-(x(j,:)-mu1)*phi;
em2=(y(j)-mu2)-(x(j,:)-mu1)*phi; 
em3=(y(j)-mu1)-(x(j,:)-mu2)*phi;
em4=(y(j)-mu2)-(x(j,:)-mu2)*phi; 
neta1=(1/sqrt(sig1))*exp(-0.5*(em1*isig1*em1'));%F(Y\S=1)
neta2=(1/sqrt(sig2))*exp(-0.5*(em2*isig2*em2'));%F(Y\S=2)
neta3=(1/sqrt(sig1))*exp(-0.5*(em3*isig3*em3'));%F(Y\S=3)
neta4=(1/sqrt(sig2))*exp(-0.5*(em4*isig4*em4'));%F(Y\S=4)
 %%%Prediction Step%%%%
 ett10=pmat*ett11;
 %%%%Update Step%%%%
ett11=ett10.*[neta1;neta2;neta3;neta4]; %joint density F(Y,S)
fit=sum(ett11);           %Marginal density F(Y)
ett11=(ett11)/fit;    %conditional density F(S\Y) the weights of the likelihood
filter(j,1:4)=ett11';      %save filter probability ett  
lik=lik+log(fit);
end
 check=-1;
 while check<0
    %backward recursion to sample from H(S[t]\S[t+1],y)
   S=zeros(T,1);
   %time T
   p1=filter(T,1);
   p2=filter(T,2);
   p3=filter(T,3);
   p4=filter(T,4);
  p=p1/(p1+p2+p3+p4);
   u=rand(1,1);
  temp=(u>=p);
  if temp==0
      S(T)=1;
  else
     p=p2/(p2+p3+p4); 
      u=rand(1,1);
  temp=(u>=p);
     if temp==0
      S(T)=2;
     else
        p=p3/(p3+p4); 
      u=rand(1,1);
  temp=(u>=p); 
  if (temp==0)
      S(T)=3;
  else
      S(T)=4;
     end
     end
  end
  
   %time t-1 to 1

 for t=T-1:-1:1
   if S(t+1)==1       
p1=pmat(1,1)*filter(t,1);
p2=pmat(1,2)*filter(t,2);
p3=pmat(1,3)*filter(t,3);
p4=pmat(1,4)*filter(t,4);

elseif S(t+1)==2
p1=pmat(2,1)*filter(t,1);
p2=pmat(2,2)*filter(t,2);
p3=pmat(2,3)*filter(t,3);
p4=pmat(2,4)*filter(t,4);
elseif S(t+1)==3
p1=pmat(3,1)*filter(t,1);
p2=pmat(3,2)*filter(t,2);
p3=pmat(3,3)*filter(t,3);
p4=pmat(3,4)*filter(t,4);
elseif S(t+1)==4
p1=pmat(4,1)*filter(t,1);
p2=pmat(4,2)*filter(t,2);
p3=pmat(4,3)*filter(t,3);
p4=pmat(4,4)*filter(t,4);

   end
   
   %sample regime numbers
   p=p1/(p1+p2+p3+p4);
   u=rand(1,1);
  temp=(u>=p);
  if temp==0
      S(t)=1;
  else
     p=p2/(p2+p3+p4); 
      u=rand(1,1);
  temp=(u>=p);
     if temp==0
      S(t)=2;
     else
         p=p3/(p3+p4);
         u=rand(1,1);
           temp=(u>=p);
     if temp==0
         S(t)=3;
     else
         S(t)=4;
     end
  end
  end 
 end
   
if sum(((S==1)+(S==3))==1)>=ncrit && sum(((S==2)+(S==4))==1)>=ncrit
    check=1;
end
 end
 %construct State variable with two regimes
 Sstar=1-((S==1)+(S==3)); %equals 0 for regime 1, 1 for regime 2
 
 %step 2 sample the transition matrix P
    %calculate the number of regime switches
    tranmat=switchg(Sstar+1,[1;2]);
     N00=tranmat(1,1); %S(t-1)=0 S(t)=0
    N01=tranmat(1,2); %S(t-1)=0 S(t)=1
    N10=tranmat(2,1); %S(t-1)=1 S(t)=0
    N11=tranmat(2,2); %S(t-1)=1 S(t)=1
    %draw from the dirichlet density
    p0=drchrnd([N00+u00;N01+u01]);
    p=p0(1,1); %p00
    p0=drchrnd([N10+u10;N11+u11]);
    q=p0(2,1); %p11
    pmat0=[p 1-q;1-p q];
    pmat=matf(pmat0,2,1);

    
    %step 3 sample AR coefficient
    mut=(S==1).*mu1+(S==2).*mu2+(S==3).*mu1+(S==4).*mu2;
    mutlag=(S==1).*mu1+(S==2).*mu1+(S==3).*mu2+(S==4).*mu2;
    sigall=(S==1).*sig1+(S==2).*sig2+(S==3).*sig1+(S==4).*sig2;
    ystar=(y-mut)./sqrt(sigall);
    xstar=(x-mutlag)./sqrt(sigall);
     V=inv(inv(Sigma0)+(xstar'*xstar));
     M=V*(inv(Sigma0)*B0+xstar'*ystar);
     phi=M+(randn(1,1)*chol(V))';
     
     %draw mu
     ystar=(y-x*phi)./sqrt(sigall);
     xstar=([(Sstar==0).*(1-phi) (Sstar==1).*(1-phi)])./repmat(sqrt(sigall),1,2);
     V=inv(inv(Sigma0M)+(xstar'*xstar));
     M=V*(inv(Sigma0M)*M0+xstar'*ystar);
     mu=M+(randn(1,2)*chol(V))';
     mu1=mu(1);
     mu2=mu(2);
     
    %regime 1
    resid=(y-mut)-(x-mutlag)*phi;
    e1=resid(Sstar==0);
    e2=resid(Sstar==1);
    
    
    %step 4 sample sigma
    
    %regime 1
   
    T1=v0+rows(e1);
    D1=d0+e1'*e1;
    %draw from IG
   z0=randn(T1,1);
    z0z0=z0'*z0;
   sig1=D1/z0z0;
   %regime 2
 
    T2=v0+rows(e2);
    D2=d0+e2'*e2;
    %draw from IG
   z0=randn(T2,1);
    z0z0=z0'*z0;
   sig2=D2/z0z0;
   
   
   
   %save and impose regime identification
   if igibbs>BURN
       chck=mu1>mu2; %mean bigger in regime 1
       if chck
           out1=[out1;([phi mu1 mu2])];
           out2=[out2;([sig1 sig2 ])];
           out3=[out3;Sstar'];
           out4=[out4;[p q]];
           count=count+1;
       end
    
   end
   igibbs=igibbs+1
end
  
figure(1)
subplot(4,2,1);
hist(out1(:,1),50);
vline(B1)
title('AR Coefficient');
axis tight

subplot(4,2,3);
hist(out1(:,2),50);
vline(MU1)
title('mean regime 1');
axis tight
subplot(4,2,4);
hist(out1(:,3),50);
title('mean regime 2');
vline(MU2)
axis tight
subplot(4,2,5);
hist(out2(:,1),50);
vline(S1)
title('\sigma_{1}');
axis tight
subplot(4,2,6);
hist(out2(:,2),50);
vline(S2)
title('\sigma_{2}');
axis tight
subplot(4,2,7);
hist(out4(:,1),50);
title('P_{00}');
vline(P0(1,1))
axis tight
subplot(4,2,8);
hist(out4(:,2),50);
title('p_{11}');
vline(P0(2,2))
axis tight

figure(2)
temp=mean(out3,1);
plot(temp,'k');
hold on
plot(strue(:,2)+strue(:,4),'c')
title('Probability of Regime 2');
legend('Estimate','True')
axis tight


