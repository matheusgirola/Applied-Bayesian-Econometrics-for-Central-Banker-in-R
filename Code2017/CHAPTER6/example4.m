clear

addpath('distributions','functions','gensys','sims_Optimization');  %change these to reflect yours!

%load data
load data


%set starting values for each parameter

Theta=[0.5 1.1 0.1 1.1 0.5 0.5 0.1 0.1 0.1];

%set bounds for each parameter
bounds=zeros(length(Theta),2);
bounds(1,:)=[0.01 5];  %sigma
bounds(2,:)=[1.01 5];  %delta
bounds(3,:)=[0.01 5];  %alpha
bounds(4,:)=[1.01 5];  %omega
bounds(5,:)=[0.01 0.999];  %rho2
bounds(6,:)=[0.01 0.999];  %rho3
bounds(7,:)=[0.01 5];  %sigma 1
bounds(8,:)=[0.01 5];  %sigma 2
bounds(9,:)=[0.01 5];  % sigma 3
options = optimset('Disp','iter','Diagnostics','on','LargeScale','off',...
    'MaxFunEvals',100000,'MaxIter',500,'TolFun',1e-05,'TolX',1e-05);

% % %simplex
 [Theta1,fval] = fminsearch(@posteriorR, Theta,options,y,bounds,1);

%%%BFGS
[fh,xh,gh,H,itct,fcount,retcodeh] = csminwel('posteriorR',Theta1,eye(length(Theta)).*0.001,[],0.0001,10000,y,bounds,1);


REPS=35000;
BURN=20000;
%Metropolis Hastings algorithm
K=0.3;
P=chol(H*K);  %choleski decomposition of inverse hessian used as starting value


bold=xh';   %starting value for DSGE parameters
out=[];
outpost=[];
naccept=0;
pold=posteriorR(bold',y,bounds,0); 
for i=1:REPS
    
    %step 1 Generate new draw from random walk
    bnew=bold+(randn(1,length(bold))*P)';
    
    %step 2 Evaluate Posterior at new draw
    pnew=posteriorR(bnew',y,bounds,0);
    
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

