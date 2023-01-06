clear;
addpath('functions');
bx1=[ 0.7 -0.1 0;  -0.1 0.7 0 ];
bx2=[ 0.1 -0.1 -3;  -0.1 0.1 -10];
sigmax=[0.5 -0.1;
        -0.1 1];

T=500;
dataout=zeros(T,2);
TAR=-1;
GAM=3;
for i=2:T
    Ystar(i)=dataout(i-1,1);
    LSTAR(i)=1./(1+exp(-GAM.*(Ystar(i)-TAR)));
     e1=1-LSTAR(i);
    e2=LSTAR(i);
    dataout(i,:)=e1.*(([dataout(i-1,:) 1]*bx1'))+...
        e2.*(([ dataout(i-1,:) 1]*bx2'))+randn(1,2)*chol(sigmax);
end


data=dataout;
L=1;
lamdaP=1;
tauP=10*lamdaP;
epsilonP=1/1000;
tarvar=1; %threshold variable
tard=1;  %delay
gammean=1; %prior mean Gamma
gamvariance=10; %prior variance Gamma
tarvariance=10; %prior variance Threshold
REPS=20000;
BURN=15000;

MaxTrys=1000;

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
Ystar=lag0(Y(:,tarvar),tard(1));


Y=Y(max([L,tard(1)])+1:end,:);
X=X(max([L,tard(1)])+1:end,:);
Ystar=Ystar(max([L,tard(1)])+1:end,:);
tarmean=mean(Ystar);  %mean of the prior on the threshold is the mean value of the threshold variable

% Additional priors for VAR coefficients
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



  Y0=[Y;yd];
  X0=[X;xd];
  
  sigma1=eye(N); %starting value for sigma
  
  beta0=vec(X0\Y0);
  beta01=beta0;
  beta02=beta0;
  beta2=beta02;
  tar=tarmean; %initial value of the threshold
  gam=gammean; %initial value of smoothness parameter
  paramold=[tar;gam]; %starting value for MH step
  paramvar=eye(2).*0.001;%diag(abs(paramold)).*2; %scaling for MH step

 
  
  naccept=0;
  igibbs=1;
  jgibbs=0;
while jgibbs<REPS-BURN
	
    
    
    %step 1: Seperate into two regimes
    %evaluate LSTAR function
    LSTAR=1./(1+exp(-gam.*(Ystar-tar)));
    e1=(1-LSTAR);
    e2=LSTAR;
    X1=X.*repmat(e1,1,cols(X));
    X2=X.*repmat(e2,1,cols(X));
    
    %step 2 Sample Coefficients and variance regime 1
    Ystarx=Y-X2*reshape(beta2,N*L+1,N);
    Y0=[Ystarx;yd];
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
    %step 3 Sample Coefficients and variance in regime 2
    Ystarx=Y-X1*reshape(beta1,N*L+1,N);
     Y0=[Ystarx;yd];
    X0=[X2;xd];
  %conditional mean of the VAR coefficients
  mstar2=vec(X0\Y0);  %ols on the appended data
  xx=X0'*X0;
  ixx2=xx\eye(cols(xx));
   [ beta2,PROBLEM2] = getcoef( mstar2,sigma1,ixx2,MaxTrys,N,L );
     if PROBLEM2
         beta2=beta02;
     else
         beta02=beta2;
     end

%     draw covariance
     e=Y0-X0*reshape(beta2,N*L+1,N);
    scale=e'*e;
    sigma1=iwpQ(rows(Y0),inv(scale));  
    
    
    
    %step 4 Sample Threshold via a Random Walk Metropolis Step
    
    paramnew=paramold+(randn(1,2)*chol(paramvar))';

     %compute conditional posterior at the old and new draw

      postnew=getvarpostx(Y,X,beta1,beta2,sigma1,L,[paramnew],tarmean,tarvariance,gammean,gamvariance,Ystar,ncrit);
      postold=getvarpostx(Y,X,beta1,beta2,sigma1,L,[paramold],tarmean,tarvariance,gammean,gamvariance,Ystar,ncrit);
 
     accept=exp(postnew-postold);
    %end
     u=rand(1,1);
     if u<accept
         paramold=paramnew;
         naccept=naccept+1;
     end
     tar=paramold(1);
     gam=paramold(2);
     arate=naccept/igibbs;
     
     if igibbs>100 && igibbs<1000
         if arate<0.2
             paramvar=paramvar*0.99999;
         elseif arate>0.5
             paramvar=paramvar*1.01;
         end
     end
     
    
    % Display progress:
   
        disp(sprintf('     Replication %s of %s acceptance %s. ', num2str(igibbs),num2str(REPS),num2str(arate))) 
   
    
     

    

     if igibbs>BURN && ~PROBLEM1 &&~ PROBLEM2
         jgibbs=jgibbs+1;
       
   bsave1(jgibbs,:)=beta1';
   bsave2(jgibbs,:)=beta2';
   sigmaS1(jgibbs,:,:)=sigma1;
  
   tsave(jgibbs,:)=[tar gam ];

     end
     igibbs=igibbs+1;
end

figure(1)
subplot(1,2,1)
hist(tsave(:,1))
vline(TAR);
title('Threshold');
subplot(1,2,2)
hist(tsave(:,2))
vline(GAM);
title('\gamma');
