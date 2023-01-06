clear
addpath('functions');
REPS=5000;
BURN=3000;
[data,names]=xlsread('\data\usdata1.xls'); %load US data
N=cols(data);

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
%priors for VAR coefficients
lamdaP=1;  %This controls the tightness of the priors on the first lag
tauP=10*lamdaP;  % this controls the tightness of the priors on sum of coefficients
epsilonP=1;  % this controls tightness of the prior on the constant
muP=mean(Y)';
sigmaP=[];
deltaP=[];
for i=1:N
    ytemp=Y(:,i);
    xtemp=[lag0(ytemp,1) ones(rows(ytemp),1)];
    ytemp=ytemp(2:end,:);
    xtemp=xtemp(2:end,:);
    btemp=xtemp\ytemp;
    etemp=ytemp-xtemp*btemp;
    stemp=etemp'*etemp/rows(ytemp);
    deltaP=[deltaP;btemp(1)];
    sigmaP=[sigmaP;stemp];
end
%dummy data to implement priors see http://ideas.repec.org/p/ecb/ecbwps/20080966.html
[yd,xd] = create_dummies(lamdaP,tauP,deltaP,epsilonP,L,muP,sigmaP,N);

%yd and xd are the dummy data. Append this to actual data
  Y0=[Y;yd];
  X0=[X;xd];
  %conditional mean of the VAR coefficients
  mstar=vec(X0\Y0);  %ols on the appended data
  xx=X0'*X0;
  ixx=xx\eye(cols(xx));  %inv(X0'X0) to be used later in the Gibbs sampling algorithm
  sigma=eye(N); %starting value for sigma
  out=zeros(REPS-BURN,36,N);
  jj=1;
  for i=1:REPS
       vstar=kron(sigma,ixx);
       beta=mstar+(randn(1,N*(N*L+1))*chol(vstar))';
       
       %draw covariance
       e=Y0-X0*reshape(beta,N*L+1,N);
    scale=e'*e;
    sigma=IWPQ(T+rows(yd),inv(scale));
    
    if i>=BURN
        %impose sign restrictions
        chck=-1;
        while chck<0
        K=randn(N,N);
        Q=getqr(K);
        A0hat=chol(sigma);
        A0hat1=(Q*A0hat);  %candidate draw
        %check signs
        e1=A0hat1(1,1)>0;  %Response of R
        e2=A0hat1(1,2)<0;  %Response of Y
        e3=A0hat1(1,3)<0;  %Response of Inflation
        e4=A0hat1(1,4)<0;  %Response of consumption
        e5=A0hat1(1,5)>0;  %Response of U
        e6=A0hat1(1,6)<0;  %Response of investment
        e7=A0hat1(1,8)<0; %response of money
        if e1+e2+e3+e4+e5+e6+e7==7
            chck=10;
        else
            %check signs but reverse them
        e1=-A0hat1(1,1)>0;  %Response of R
        e2=-A0hat1(1,2)<0;  %Response of Y
        e3=-A0hat1(1,3)<0;  %Response of Inflation 
        e4=-A0hat1(1,4)<0;  %Response of consumption
        e5=-A0hat1(1,5)>0;  %Response of U
        e6=-A0hat1(1,6)<0;  %Response of investment
        e7=-A0hat1(1,8)<0; %response of money
        if e1+e2+e3+e4+e5+e6+e7==7
            A0hat1(1,1:N)=-A0hat1(1,1:N);
            chck=10;
        end
        end
        end
        
yhat=zeros(36,N);
vhat=zeros(36,N);
vhat(3,1)=1; %shock to the Federal Funds rate

for j=3:36
 yhat(j,:)=[ yhat(j-1,:) yhat(j-2,:) 0]*reshape(beta,N*L+1,N)+vhat(j,:)*A0hat1;
end

 out(jj,:,:)=yhat;
 jj=jj+1;
    end
   
  end
  
  figure(1);
  subplot(4,3,1)
  temp=out(:,:,1);
  temp1=squeeze(prctile(temp,[50 16 84],1))';
  plot(temp1(3:end,:));
  title('FederalFunds Rate');
  axis tight
  
  
  
  subplot(4,3,2)
  temp=out(:,:,2);
  temp1=squeeze(prctile(temp,[50 16 84],1))';
  plot(temp1(3:end,:));
  title('GDP Growth');
  axis tight
  
  
  subplot(4,3,3)
  temp=out(:,:,3);
  temp1=squeeze(prctile(temp,[50 16 84],1))';
  plot(temp1(3:end,:));
  title('CPI Inflation');
  axis tight
  
  
   subplot(4,3,4)
  temp=out(:,:,4);
  temp1=squeeze(prctile(temp,[50 16 84],1))';
  plot(temp1(3:end,:));
  title('PCE Growth');
  axis tight
  
  subplot(4,3,5)
  temp=out(:,:,5);
  temp1=squeeze(prctile(temp,[50 16 84],1))';
  plot(temp1(3:end,:));
  title('Unemployment');
  axis tight
  
   subplot(4,3,6)
  temp=out(:,:,6);
  temp1=squeeze(prctile(temp,[50 16 84],1))';
  plot(temp1(3:end,:));
  title('Investment');
  axis tight
  
  
   subplot(4,3,7)
  temp=out(:,:,7);
  temp1=squeeze(prctile(temp,[50 16 84],1))';
  plot(temp1(3:end,:));
  title('Net Exports');
  axis tight
  
  
   subplot(4,3,8)
  temp=out(:,:,8);
  temp1=squeeze(prctile(temp,[50 16 84],1))';
  plot(temp1(3:end,:));
  title('M2');
  axis tight
  
  
   subplot(4,3,9)
  temp=out(:,:,9);
  temp1=squeeze(prctile(temp,[50 16 84],1))';
  plot(temp1(3:end,:));
  title('10 year Government Bond Yield');
  axis tight
  
  
   subplot(4,3,10)
  temp=out(:,:,10);
  temp1=squeeze(prctile(temp,[50 16 84],1))';
  plot(temp1(3:end,:));
  title('Stock Price Growth');
  axis tight
  
  
  subplot(4,3,11)
  temp=out(:,:,11);
  temp1=squeeze(prctile(temp,[50 16 84],1))';
  plot(temp1(3:end,:));
  title('Yen Dollar Rate');
  axis tight
  
  
  
  
  
  
  
  
      