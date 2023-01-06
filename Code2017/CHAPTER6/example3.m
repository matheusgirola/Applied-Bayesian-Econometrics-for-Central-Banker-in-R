clear

addpath('distributions','functions','gensys','sims_Optimization');  %change these to reflect yours!

%load data
load data


%set starting values for each parameter

Theta=[0.5 1.1 0.1 1.1 0.5 0.5 0.5 0.1 0.1 0.1];

%set bounds for each parameter
bounds=zeros(length(Theta),2);
bounds(1,:)=[0.01 5];  %sigma
bounds(2,:)=[1.01 5];  %delta
bounds(3,:)=[0.01 5];  %alpha
bounds(4,:)=[1.01 5];  %omega
bounds(5,:)=[0.01 0.999]; %rho1
bounds(6,:)=[0.01 0.999];  %rho2
bounds(7,:)=[0.01 0.999];  %rho3
bounds(8,:)=[0.01 5];  %sigma 1
bounds(9,:)=[0.01 5];  %sigma 2
bounds(10,:)=[0.01 5];  % sigma 3
options = optimset('Disp','iter','Diagnostics','on','LargeScale','off',...
    'MaxFunEvals',100000,'MaxIter',500,'TolFun',1e-05,'TolX',1e-05);

% % %simplex
 [Theta1,fval] = fminsearch(@posterior, Theta,options,y,bounds,1);

%%%BFGS
[fh,xh,gh,H,itct,fcount,retcodeh] = csminwel('posterior',Theta1,eye(length(Theta)).*0.001,[],0.0001,10000,y,bounds,1);


REPS=35000;
BURN=20000;
%Metropolis Hastings algorithm
K=0.3;
P=chol(H*K);  %choleski decomposition of inverse hessian used as starting value


bold=xh';   %starting value for DSGE parameters
out=[];
outpost=[];
naccept=0;
pold=posterior(bold',y,bounds,0); 
for i=1:REPS
    
    %step 1 Generate new draw from random walk
    bnew=bold+(randn(1,length(bold))*P)';
    
    %step 2 Evaluate Posterior at new draw
    pnew=posterior(bnew',y,bounds,0);
    
%     %step 3 Evaluate Posterior at old draw
     
    %compute accpetance probability
    
    if pnew==-inf;
        accept=0;
    else
        accept=min([exp(pnew-pold) 1]);
    end
    
    u=rand(1,1);
    if u<accept
        bold=bnew;
        pold=pnew;
        naccept=naccept+1;
    end
    if i>BURN
out=[out;bold'];
outpost=[outpost;pold];
    end
    arate=naccept/i;
        disp(sprintf(' Replication %s of %s., acceptance %s', ... 
             num2str(i ), num2str(REPS),num2str(arate)) );
end

draws=out;


'----------------------Posterior Mean--------------------------------'
mean(draws)
'----------------------Standard Deviation-----------------------------'
std(draws)
'----------------------Lower Bound-------------------------------------'
prctile(draws,16)
'----------------------Upper Bound-------------------------------------'
prctile(draws,84)
'-----------------------------------------------------------------------'
'Acceptance Rate';naccept/REPS;

outp=simprior(size(draws,1));















figure(1)
subplot(5,2,1);
[xx,ff]=ksdensity(draws(:,1));
[xx1,ff1]=ksdensity(outp(:,1));
plot(ff,xx,ff1,xx1)
vline(Thetat(:,1))
vline(xh(1),'c')
title('\sigma');

subplot(5,2,2);
[xx,ff]=ksdensity(draws(:,2));
[xx1,ff1]=ksdensity(1+outp(:,2));
plot(ff,xx,ff1,xx1)
vline(Thetat(:,2))
vline(xh(2),'c')
title('\delta');


subplot(5,2,3);
[xx,ff]=ksdensity(draws(:,3));
[xx1,ff1]=ksdensity(outp(:,3));
plot(ff,xx,ff1,xx1)
vline(Thetat(:,3))
vline(xh(3),'c')
title('\alpha');


subplot(5,2,4);
[xx,ff]=ksdensity(draws(:,4));
[xx1,ff1]=ksdensity(1+outp(:,4));
plot(ff,xx,ff1,xx1)
vline(Thetat(:,4))
vline(xh(4),'c')
title('\omega');


subplot(5,2,5);
[xx,ff]=ksdensity(draws(:,5));
[xx1,ff1]=ksdensity(outp(:,5));
plot(ff,xx,ff1,xx1)
vline(Thetat(:,5))
vline(xh(5),'c')
title('\rho_1');


subplot(5,2,6);
[xx,ff]=ksdensity(draws(:,6));
[xx1,ff1]=ksdensity(outp(:,6));
plot(ff,xx,ff1,xx1)
vline(Thetat(:,6))
vline(xh(6),'c')
title('\rho_2');


subplot(5,2,7);
[xx,ff]=ksdensity(draws(:,7));
[xx1,ff1]=ksdensity(outp(:,7));
plot(ff,xx,ff1,xx1)
vline(Thetat(:,7))
vline(xh(7),'c')
title('\rho_3');


subplot(5,2,8);
[xx,ff]=ksdensity(draws(:,8));
[xx1,ff1]=ksdensity(outp(:,8));
plot(ff,xx,ff1,xx1)
vline(Thetat(:,8))
vline(xh(8),'c')
title('\sigma_1');


subplot(5,2,9);
[xx,ff]=ksdensity(draws(:,9));
[xx1,ff1]=ksdensity(outp(:,9));
plot(ff,xx,ff1,xx1)
vline(Thetat(:,9))
vline(xh(9),'c')
title('\sigma_2');

subplot(5,2,10);
[xx,ff]=ksdensity(draws(:,10));
[xx1,ff1]=ksdensity(outp(:,10));
plot(ff,xx,ff1,xx1)
vline(Thetat(:,10))
vline(xh(10),'c')
title('\sigma_3');
legend('Posterior','Prior');



%compute impulse responses to a policy shock

outy=[];
outp=[];
outi=[];
for ii=1:size(draws,1)
    [  PP, QQ, PROBLEM ] = model_solve( draws(ii,:) );

% Generate data 
Sigma=diag(draws(ii,8:10));
t=20;
ee=zeros(t,3);
ee(2,3)=1;
ee=ee*chol(Sigma);
ehat=zeros(t,6);
ymat=zeros(t,6);
for i=2:t
    %iid shocks
    ehat(i,:)=(QQ*ee(i,:)')';
    %data
    ymat(i,:)=(PP*ymat(i-1,:)')'+[ehat(i,:)];
end
outy=[outy ymat(:,1)];
outp=[outp ymat(:,2)];
outi=[outi ymat(:,3)];
end
load IRFT %true IRF
figure(2);
subplot(2,2,1);
tmp=prctile(outy,[50 16 84],2);
plot(tmp,'r');
hold on
plot(IRFT(:,1),'.-k')
axis tight
title('Y');
subplot(2,2,2);
tmp=prctile(outp,[50 16 84],2);
plot(tmp,'r');
hold on
plot(IRFT(:,2),'.-k')
axis tight
title('P');

subplot(2,2,3);
tmp=prctile(outi,[50 16 84],2);
plot(tmp,'r');
hold on
plot(IRFT(:,3),'.-k')
axis tight
title('R');

%calculate the marginal likelihood using Gelfand and Dey method
%posterior mean and variance
pmean=mean(draws);
pvar=cov(draws);
ipvar=inv(pvar);
dpvar=log(det(pvar));
lpost_mode=max(outpost);
p=0.1; %critical value of the Chi-squared distribution
npara=size(draws,2); %number of parameters
critval = chi2inv(p,npara);
    tmp = 0;
    for i = 1:size(draws,1);
        deviation  = (draws(i,:)-pmean)*ipvar*((draws(i,:)-pmean))';
        if deviation <= critval;
            lftheta = -log(p)-(npara*log(2*pi)+dpvar+deviation)/2;
            tmp = tmp + exp(lftheta - outpost(i)+lpost_mode);
        end;    
    end;
mlik=lpost_mode-log(tmp/size(draws,1));

disp('Gelfand and Dey log Marginal Likelihood');
disp(num2str(mlik));

