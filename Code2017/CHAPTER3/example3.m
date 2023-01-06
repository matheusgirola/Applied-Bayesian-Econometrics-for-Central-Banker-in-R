clear
addpath('functions');


% a TVP-VAR with dlog(GDP) dlog(CPI) and R for the US 1962 2004
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

%initialise 
Q=Q0;
R=sigma0;


%remove intial Sample
Y=Y(T0+1:end,:);
X=X(T0+1:end,:);
T=rows(X);

%Gibbs sampling algorithm Step 2
reps=110000;
burn=109000;
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

%step4 sample R from the IW distribution
scaleR=(error'*error);
R=iwpQ(T,inv(scaleR));
if m>burn
    %save output from Gibbs sampler
    out1(mm,1:T,:)=beta2;
    out2(mm,1:N,1:N)=R;
    out3(mm,1:N*(N*L+1),1:N*(N*L+1))=Q;
    mm=mm+1;
end

end
%save results
save tvp.mat out1 out2 out3


%compute irf to a policy shock using sign restrictions
horz=40;% impulse response horizon
irfmat=zeros(size(out1,1),T,horz,N); %empty matrix to save impulse response to a policy shock
for i=1:size(out1,1);
    sigma=squeeze(out2(i,:,:));
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
    for j=1:size(out1,2)
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
subplot(2,2,1);
mesh(TT,HH,irf1')
ylabel('Impulse Horizon');
xlabel('Time');
axis tight
title('GDP growth');

subplot(2,2,2);
mesh(TT,HH,irf2')
ylabel('Impulse Horizon');
xlabel('Time');
axis tight
title('Inflation');

subplot(2,2,3);
mesh(TT,HH,irf3')
ylabel('Impulse Horizon');
xlabel('Time');
axis tight
title('Federal Funds Rate');






