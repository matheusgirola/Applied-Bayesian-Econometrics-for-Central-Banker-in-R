clear
addpath('functions');


% a TVP-VAR with stochastic volatility using dlog(GDP) dlog(CPI) and R for the US 1962 2004
%load data
data=xlsread('\data\usdata.xls')/100;

N=size(data,2);
L=2;   %number of lags in the VAR
Y=data;
X=[ lag0(Y,1) lag0(Y,2) ones(size(Y,1),1) ];
Y=Y(3:end,:);
X=X(3:end,:);


%step 1 set starting values and priors using a pre-sample of 10 years
T0=40;
y0=Y(1:T0,:);
x0=X(1:T0,:);
b0=x0\y0;
e0=y0-x0*b0;
sigma0=(e0'*e0)/T0;

V0=kron(sigma0,inv(x0'*x0));


%priors for the variance of the transition equation
Q0=V0*T0*3.5e-04;  %prior for the variance of the transition equation error
P00=V0;                    % variance of the intial state vector  variance of state variable p[t-1/t-1]
beta0=vec(b0)';       % intial state vector   %state variable  b[t-1/t-1]

%priors and starting values for aij
C0=chol(sigma0);
C0=C0./repmat(diag(C0),1,N);
C0=inv(C0)';

a10=C0(2,1);   %intial state vector
a20=C0(3,1:2);  %intial state vector second equation
pa10=abs(a10)*10;  %variance of the state vector
pa20=diag(a20)*10;  %variance of the state vector
D10=10^(-3);  %prior scale matrix for D1
D20=10^(-3)*eye(2); %prior scale matrix for D2



%remove intial Sample
Y=Y(T0+1:end,:);
X=X(T0+1:end,:);
T=rows(X);


%priors and starting values for the stochastic vol
hlast=(diff(Y).^2)+0.0001;
hlast=[hlast(1:2,:);hlast];  %rough intial guess for svol
g=ones(3,1);  %rough guess for the variance of the transition equation
g0=0.01^2;  %scale parameter for inverse gamma
Tg0=1;
mubar=log(diag(sigma0));
sigmabar=10;

%initialise parameters
Q=Q0;
D1=D10;
D2=D20;
a1=repmat(a10,T,1);
a2=repmat(a20,T,1);


%Gibbs sampling algorithm Step 2
reps=100000;
burn=99000;
mm=1;
for m=1:reps
m
%%Step 2a Set up matrices for the Kalman Filter

ns=cols(beta0);
F=eye(ns);
mu=0;
beta_tt=[];          %will hold the filtered state variable
ptt=zeros(T,ns,ns);    % will hold its variance
beta11=beta0;
p11=P00;


% %%%%%%%%%%%Step 2b run Kalman Filter

for i=1:T
   x=kron(eye(N),X(i,:));

a=[a1(i) a2(i,:)];
A=chofac(N,a');
H=diag(hlast(i+1,:));
R=inv(A)*H*inv(A)';
    %Prediction
beta10=mu+beta11*F';
p10=F*p11*F'+Q;
yhat=(x*(beta10)')';                                               
eta=Y(i,:)-yhat;
feta=(x*p10*x')+R;
%updating
K=(p10*x')*inv(feta);
beta11=(beta10'+K*eta')';
p11=p10-K*(x*p10);

ptt(i,:,:)=p11;
beta_tt=[beta_tt;beta11];

end


%%%%%%%%%%%end of Kalman Filter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%step 2c Backward recursion to calculate the mean and variance of the distribution of the state
%vector
chck=-1;
while chck<0
beta2 = zeros(T,ns);   %this will hold the draw of the state variable
wa=randn(T,ns);
error=zeros(T,N);
roots=zeros(T,1);

i=T;  %period t
p00=squeeze(ptt(i,:,:)); 
beta2(i,:)=beta_tt(i:i,:)+(wa(i:i,:)*chol(p00));   %draw for beta in period t from N(beta_tt,ptt)
error(i,:)=Y(i,:)-X(i,:)*reshape(beta2(i:i,:),N*L+1,N);  %var residuals
roots(i)=stability(beta2(i,:)',N,L);

%periods t-1..to .1
for i=T-1:-1:1
pt=squeeze(ptt(i,:,:));
bm=beta_tt(i:i,:)+(pt*F'*inv(F*pt*F'+Q)*(beta2(i+1:i+1,:)-beta_tt(i,:)*F')')';  %update the filtered beta for information contained in beta[t+1]                                                                                 %i.e. beta2(i+1:i+1,:) eq 8.16 pp193 in Kim Nelson
pm=pt-pt*F'*inv(F*pt*F'+Q)*F*pt;  %update covariance of beta
beta2(i:i,:)=bm+(wa(i:i,:)*chol(pm));  %draw for beta in period t from N(bm,pm)eq 8.17 pp193 in Kim Nelson
error(i,:)=Y(i,:)-X(i,:)*reshape(beta2(i:i,:),N*L+1,N);  %var residuals
roots(i)=stability(beta2(i,:)',N,L);
end

if sum(roots)==0
    chck=1;
end
end

% step 3 sample Q from the IW distribution
errorq=diff(beta2);
scaleQ=(errorq'*errorq)+Q0;
Q=iwpQ(T+T0,inv(scaleQ));

%step4 sample aij using the carter kohn algorithm
v3=error(:,3);
v2=error(:,2);
v1=error(:,1);

[a1,trash]=carterkohn1(a10,pa10,hlast(:,2),D1,v2,-v1);
[a2,trash]=carterkohn1(a20,pa20,hlast(:,3),D2,v3,[-v1 -v2]);

%step 5 sample D1 and D2
a1errors=diff(a1);
D1=IG(T0,D10,a1errors);  %draw from the inverse Gamma distribution
a2errors=diff(a2);
scaleD2=(a2errors'*a2errors)+D20;
D2=iwpQ(T+T0,inv(scaleD2)); %draw from inverse Wishart


%step 6 sample h_i seperately for i=1,3
 %step 6a calculate epsilon=A*v
 epsilon=[];
 for i=1:T
     a=[a1(i) a2(i,:)];
     A=chofac(N,a');
     epsilon=[epsilon;error(i,:)*A'];
 end
 %sample stochastic vol for each epsilon using the MH algorithm
  hnew=[];
  for i=1:N
      htemp=getsvol(hlast(:,i),g(i),mubar(i),sigmabar,epsilon(:,i));
      hnew=[hnew htemp];
  end
hlast=hnew;

%step 7 Sample G for IG distribution
for i=1:N
    gerrors=diff(log(hnew(:,i)));
g(i)=IG(Tg0,g0,gerrors);  %draw from the inverse Gamma distribution
end




if m>burn
    %save output from Gibbs sampler
    out1(mm,1:T,:)=beta2;
    out2(mm,1:T,1:N)=hlast(2:end,:);
    out3(mm,1:N*(N*L+1),1:N*(N*L+1))=Q;
    out4(mm,1:T,1:(N*(N-1))/2)=[a1 a2];
    out5(mm,1)=D1;
    out6(mm,1:2,1:2)=D2;
    out7(mm,1:N)=g';
    mm=mm+1;
end

end
%save results
save tvp.mat out1 out2 out3 out4 out5 out6 out7


%compute irf to a policy shock using sign restrictions
horz=40;% impulse response horizon
irfmat=zeros(size(out1,1),T,horz,N); %empty matrix to save impulse response to a policy shock
for i=1:size(out1,1);
    
    for j=1:size(out1,2)
       
       H=diag(squeeze(out2(i,j,:)));
       a=squeeze(out4(i,j,:));
        A=chofac(N,a);
        sigma=inv(A)*H*inv(A)';  %covariance matrix
   %sign restrictions
        chck=-1;
        while chck<0
        K=randn(N,N);
        QQ=getQR(K);
        A0hat=chol(sigma);
        A0hat1=(QQ*A0hat);  %candidate draw
        for m=1:N
        %check signs in each row
        e1=A0hat1(m,1)<0;  %Response of Y
        e2=A0hat1(m,2)<0;  %Response of P
        e3=A0hat1(m,3)>0;  %Response of R
        
        if e1+e2+e3==3
            MP=A0hat1(m,:);
            chck=10;
        else
            %check signs but reverse them
         e1=-A0hat1(m,1)<0;  %Response of Y
        e2=-A0hat1(m,2)<0;  %Response of P
        e3=-A0hat1(m,3)>0;  %Response of R
        
        if e1+e2+e3==3
            MP=-A0hat1(m,:);
            chck=10;
        end
        end
        end
        end
        %re-shuffle rows of A0hat1 and insert MP in the first row
        A0x=[]; %will hold rows of A0hat1 not equal to MP
        for m=1:N
            ee=sum(abs(A0hat1(m,:))==abs(MP));
            if ee==0
                A0x=[A0x;A0hat1(m,:)];
            end
        end
        A0new=[A0x;MP]; %A0 matrix to be used for impulse response
    shock=[0 0 1]; 
        btemp=squeeze(out1(i,j,:));
        btemp=reshape(btemp,N*L+1,N);
        zz=irfsim(btemp,N,L,A0new,shock,horz+L);
        zz=zz./repmat(zz(1,3),horz,N);
        irfmat(i,j,:,:)=zz;
    end
end
TT=1964.75:0.25:2010.5;
HH=0:horz-1;
irf1=squeeze(median(irfmat(:,:,:,1),1));
irf2=squeeze(median(irfmat(:,:,:,2),1));
irf3=squeeze(median(irfmat(:,:,:,3),1));

figure(1)
subplot(2,3,1);
mesh(TT,HH,irf1')
ylabel('Impulse Horizon');
xlabel('Time');
axis tight
title('GDP growth');

subplot(2,3,2);
mesh(TT,HH,irf2')
ylabel('Impulse Horizon');
xlabel('Time');
axis tight
title('Inflation');

subplot(2,3,3);
mesh(TT,HH,irf3')
ylabel('Impulse Horizon');
xlabel('Time');
axis tight
title('Federal Funds Rate');

temp=(squeeze(out2(:,:,1)));
subplot(2,3,4);
plot(TT,prctile(temp,[50 16 84],1));
axis tight
title('Stochastic Volatility GDP Growth')


temp=(squeeze(out2(:,:,2)));
subplot(2,3,5);
plot(TT,prctile(temp,[50 16 84],1));
axis tight
title('Stochastic Volatility Inflation')


temp=(squeeze(out2(:,:,3)));
subplot(2,3,6);
plot(TT,prctile(temp,[50 16 84],1));
axis tight
title('Stochastic Volatility Federal Funds Rate')



%recursive means
mbeta=rmean1(out1,20);
mh=rmean1(out2,20);
ma=rmean1(out4,20);

figure(2)
subplot(2,2,1);
mesh(mbeta);
axis tight
title('\beta_{t}');
xlabel('Vectorised parameters');
ylabel('Draws')

subplot(2,2,2);
mesh(mh);
axis tight
title('h_{it}');
xlabel('Vectorised parameters');
ylabel('Draws')


subplot(2,2,3);
mesh(ma);
axis tight
title('a_{ijt}');
xlabel('Vectorised parameters');
ylabel('Draws')





