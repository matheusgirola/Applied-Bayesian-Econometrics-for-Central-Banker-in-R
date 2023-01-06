function[beta2,error,roots,problem]=carterkohnvar(Y,X,Q,iamat,hlast,beta0,P00,L,CHECK,maxdraws,EX)
                                     %carterkohnvarold(Y,X,Q,iamat,hlast,beta0,P00,L,CHECK,maxdraws,EX)
T=rows(Y);
N=cols(Y);
%%Step 2a Set up matrices for the Kalman Filter

ns=cols(beta0);
F=eye(ns);
mu=0;
beta11=beta0;
p11=P00;


% %%%%%%%%%%%Step 2b run Kalman Filter

for i=1:T
   x=kron(eye(N),X(i,:));
H=diag(hlast(i+1,:));
R=iamat*H*iamat';
    %Prediction
beta10=beta11;
p10=p11;
yhat=(x*(beta10)')';                                               
eta=Y(i,:)-yhat;
feta=(x*p10*x')+R;
%updating
K=(p10*x')*invpd(feta);
beta11=(beta10'+K*eta')';
p11=p10-K*(x*p10);

ptt=p11;
beta_tt=beta11;

end


%%%%%%%%%%%end of Kalman Filter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%step 2c Backward recursion to calculate the mean and variance of the distribution of the state
%vector
chck=-1;
problem=0;
trys=1;
while chck<0 && trys<=maxdraws
    
wa=randn(1,ns);


beta2=beta_tt+(wa(1,:)*cholx(ptt));   %draw for beta in period t from N(beta_tt,ptt)
error=Y-X*reshape(beta2,N*L+EX,N);  %var residuals

roots=stability(beta2',N,L,EX);


if CHECK
if sum(roots)==0
    chck=1;
else
    trys=trys+1;
end
else
 chck=1;
end
end
if CHECK
    if chck<0
        problem=1;
    end
end