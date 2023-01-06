clear
addpath('sims_Optimization');
addpath('functions');


%create artificial data for time-varying paramteer model
T=300;
%generate artificial data on a time-varying parameter model
N=2;
Q=eye(N,N)*0.1;
R=2;
F(1,1)=0.95;
F(2,2)=0.95;
e=randn(T,2);
v=randn(T,1);
x=[randn(T,1) ones(T,1)];
y=zeros(T,1);
b=zeros(T,2);
MU=[0.1 -0.1];

for j=2:T
    b(j,:)=MU+b(j-1,:)*F'+e(j,:)*chol(Q);
    y(j,:)=x(j,:)*b(j,:)'+v(j,:)*sqrt(R);
end
  
TRUE=[diag(F);MU';R;diag(Q)];
%***********step 1 set priors for each parameter
% F~N(F0,VF0)
F0=ones(2,1);
VF0=eye(2)*0.2;
%MU~N(MU0,VMU0);
MU0=zeros(2,1);
VMU0=eye(2);
%1/R~Gamma(R0,VR0)
R0=1;
VR0=1;
%1/Q(i,i)~Gamma(Q0,VQ0)
Q0=0.1;
VQ0=1;


%***************step 2 estimate model via maximum likelihood
theta0=ones(7,1).*0.1;
[FF,AA,gh,hess,itct,fcount,retcodeh] = csminwel('posterior',theta0,eye(length(theta0))*.1,[],1e-15,1000,y,x,F0,VF0,MU0,VMU0,R0,VR0,Q0,VQ0);

%**************step 2 set scale factor for the metropolis hastings
K=0.4;  %scaling factor
P=(chol(hess*K)); %compute variance of the random walk


Gammaold=AA;
REPS=5000;
BURN=3000;
naccept=0;
out1=zeros(REPS-BURN,7);
out2=zeros(REPS-BURN,1);

%compute posterior at old draw
           %compute -1*likelihood at old draw
        lik=likelihoodTVP(Gammaold,y,x);
        %evaluate prior for each set of parameters
        F=Gammaold(1:2);
        MU=Gammaold(3:4);
        R=Gammaold(5);
        Q=Gammaold(6:7);
        
        Fprior=log(mvnpdf(F,F0,VF0));
        %prior for MU
        MUprior=log(mvnpdf(MU,MU0,VMU0));
        %prior for 1/R
        Rprior=gampdf1(VR0,R0,1/R);
        %prior for 1/Q
         Qprior=0;
        for i=1:2
         Qprior=Qprior+(gampdf1(VQ0,Q0,1/Q(i))); 
        end
        %joint prior is the sum of these
        priorold=Fprior+MUprior+Rprior+Qprior;
        posteriorOLD=-lik+priorold;
        jj=1;
for j=1:REPS

    %step 1 draw new Gamma
    Gammanew=Gammaold+(randn(1,7)*P)';
    
    %step 2 check elements of D are positive
    check=sum(Gammanew(5:end)<0) && sum(Gammanew(1:2)>1);
    if check
         posteriorNEW=-1000000;
    else
        %compute -1*likelihood at new draw
        lik=likelihoodTVP(Gammanew,y,x);
       F=Gammanew(1:2);
        MU=Gammanew(3:4);
        R=Gammanew(5);
        Q=Gammanew(6:7);
        
        Fprior=log(mvnpdf(F,F0,VF0));
        %prior for MU
        MUprior=log(mvnpdf(MU,MU0,VMU0));
        %prior for 1/R
        Rprior=gampdf1(VR0,R0,1/R);
        %prior for 1/Q
         Qprior=0;
        for i=1:2
         Qprior=Qprior+(gampdf1(VQ0,Q0,1/Q(i))); 
        end
        %joint prior is the sum of these
        
        priornew=Fprior+MUprior+Rprior+Qprior;
        posteriorNEW=-lik+priornew;
    end
        
        
        accept=min([exp(posteriorNEW-posteriorOLD);1]);   %min(accept,1)
    
    u=rand(1,1);  %random number from the uniform dist
    
    if u<accept
        Gammaold=Gammanew;  %accept draw
        posteriorOLD=posteriorNEW;
        naccept=naccept+1;  %count number of acceptances  
    end
      
     
      ARATE=naccept/j;
      if j>500 && j<1500
      if ARATE > 0.4;
          P=P*1.00000001;
      elseif ARATE<0.21;
          P=P*0.99;
      end
      end
      if j>BURN
      out1(jj,:)=Gammaold';
      out2(jj,:)=posteriorOLD;
      jj=jj+1;
      end
end


subplot(3,3,1);
plot([out1(:,1) repmat(TRUE(1),size(out1,1),1)]);
title('F_{1}');
subplot(3,3,2);
plot([out1(:,2) repmat(TRUE(2),size(out1,1),1)]);
title('F_{2}');
subplot(3,3,3);
plot([out1(:,3) repmat(TRUE(3),size(out1,1),1)]);
title('\mu_{1}');
subplot(3,3,4);
plot([out1(:,4) repmat(TRUE(4),size(out1,1),1)]);
title('\mu_{2}');
subplot(3,3,5);
plot([out1(:,5) repmat(TRUE(5),size(out1,1),1)]);
title('R');
subplot(3,3,6);
plot([out1(:,6) repmat(TRUE(6),size(out1,1),1)]);
title('Q_{1}');
subplot(3,3,7);
plot([out1(:,7) repmat(TRUE(7),size(out1,1),1)]);
title('Q_{2}');
legend('MH draws','True value');
