clear
addpath('functions');
data=xlsread('\data\tvardata.xlsx');
L=2;  %lag length
tard=2;  %delay
tarvar=3;  % threshold variable is the column number tarvar in data
MaxTrys=1000;
tarscale=0.1;  %scaling parameter for RW Metropolis algorithm
REPS=10000;
BURN=8000;
HORZ=40;  %impulse response horizon

%prepare data


Y=data;
N=cols(Y);
ncrit=(N*L+1);

%take lags
X=[];
for j=1:L
X=[X lag0(data,j) ];
end
X=[X ones(rows(X),1)];

%compute threshold variable
Ystar=lag0(Y(:,tarvar),tard);

Y=Y(max([L,tard(1)])+1:end,:);
X=X(max([L,tard(1)])+1:end,:);
Ystar=Ystar(max([L,tard(1)])+1:end,:);
tarmean=mean(Ystar);  %mean of the prior on the threshold is the mean value of the threshold variable
tarvariance=10;   %variance of the prior

% Additional priors for VAR coefficients
lamdaP  = 1;
tauP    = 10*lamdaP;
epsilonP= 1/10000;
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

T=rows(Y);
%append
  Y0=[Y;yd];
  X0=[X;xd];
  
  sigma1=eye(N); %starting value for sigma
  sigma2=eye(N);
  beta0=vec(X0\Y0);
  beta01=beta0;
  beta02=beta0;
  tar=tarmean; %initial value of the threshold 
  tarold=tar;
  naccept=0;
  
  
  
  
  %gibbs algorithm
  jgibbs=1;
  for igibbs=1:REPS
  
  
  
  
  
  %step 1: Seperate into two regimes
    e1=Ystar<=tar;
    e2=Ystar>tar;
    
    Y1=Y(e1,:);
    X1=X(e1,:);
    
    Y2=Y(e2,:);
    X2=X(e2,:);
    
    %step 2 Sample Coefficients and variance regime 1
    
    Y0=[Y1;yd];
    X0=[X1;xd];
  %conditional mean of the VAR coefficients
  mstar1=vec(X0\Y0);  %ols on the appended data
  xx=X0'*X0;
  ixx1=xx\eye(cols(xx));
   [ beta1,PROBLEM1] = getcoef( mstar1,sigma1,ixx1,MaxTrys,N,L );
     if PROBLEM1
         beta1=beta01;
     else
         beta01=beta1;
     end
     
     %draw covariance
     e=Y0-X0*reshape(beta1,N*L+1,N);
    scale=e'*e;
    sigma1=iwpQ(rows(Y0),inv(scale));  
    
    
    %step 3 Sample Coefficients and variance in regime 2
    
     Y0=[Y2;yd];
    X0=[X2;xd];
  %conditional mean of the VAR coefficients
  mstar2=vec(X0\Y0);  %ols on the appended data
  xx=X0'*X0;
  ixx2=xx\eye(cols(xx));
   [ beta2,PROBLEM2] = getcoef( mstar2,sigma2,ixx2,MaxTrys,N,L );
     if PROBLEM2
         beta2=beta02;
     else
         beta02=beta2;
     end
     
     %draw covariance
     e=Y0-X0*reshape(beta2,N*L+1,N);
    scale=e'*e;
    sigma2=iwpQ(rows(Y0),inv(scale)); 
    
    
    
    %step 4 Sample Threshold via a Random Walk Metropolis Step
    
    tarnew=tarold+randn(1,1)*sqrt(tarscale);
   

      postnew=getvarpost(Y,X,beta1,beta2,sigma1,sigma2,L,tarnew,tarmean,tarvariance,Ystar,ncrit);
      postold=getvarpost(Y,X,beta1,beta2,sigma1,sigma2,L,tarold,tarmean,tarvariance,Ystar,ncrit);

     accept=exp(postnew-postold);
     u=rand(1,1);
     if u<accept
         tarold=tarnew;
         naccept=naccept+1;
     end
     tar=tarold;
     arate=naccept/igibbs;
  if igibbs>100 && igibbs<1100
         if arate<0.2
             tarscale=tarscale*0.99;
         elseif arate>0.4
             tarscale=tarscale*1.01;
         end
  end
  
  disp(sprintf(' Replication %s of %s acceptance %s. ', ... 
             num2str(igibbs), num2str(REPS),num2str(arate) ));
         
         if igibbs>BURN
             %impulse response analysis
             A01=chol(sigma1);
             A02=chol(sigma2);
             irf1=irfsim(reshape(beta1,N*L+1,N),N,L,A01,[0 0 1],HORZ);
             irf2=irfsim(reshape(beta2,N*L+1,N),N,L,A02,[0 0 1],HORZ);
             irf1=irf1./irf1(1,3);
             irf2=irf2./irf2(1,3);
             
             %save results
             irfmat1(jgibbs,:,:)=irf1;
             irfmat2(jgibbs,:,:)=irf2;
             smat(:,jgibbs)=e1;
             jgibbs=jgibbs+1;
         end
             
             
             
  
  end
  
  
  figure(1)
  TT=1948:0.25:2012;
  subplot(4,2,[1:2])
  plot(TT,[mean(smat==1,2) Ystar])
  axis tight
  legend('S_{t}','Threshold Variable')
  
  subplot(4,2,3);
  plot(prctile(irfmat1(:,:,1),[50 16 84],1)')
  title('Regime 1');
  ylabel('GDP Growth');
  axis tight
  subplot(4,2,4);
  plot(prctile(irfmat2(:,:,1),[50 16 84],1)')
  title('Regime 2');
  
  axis tight
  
  
   subplot(4,2,5);
  plot(prctile(irfmat1(:,:,2),[50 16 84],1)')
  title('Regime 1');
  ylabel('Inflation');
  axis tight
  subplot(4,2,6);
  plot(prctile(irfmat2(:,:,2),[50 16 84],1)')
  title('Regime 2');
  
  axis tight
  
  
  subplot(4,2,7);
  plot(prctile(irfmat1(:,:,3),[50 16 84],1)')
  title('Regime 1');
  ylabel('Interest rate');
  axis tight
  subplot(4,2,8);
  plot(prctile(irfmat2(:,:,3),[50 16 84],1)')
  title('Regime 2');
  
  axis tight
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  