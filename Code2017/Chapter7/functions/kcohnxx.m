function [beta2,error,beta_tt]=kcohnxx(Y,X,a1,hlast,beta0,P00,L)
Q=0;
%%Step 2a Set up matrices for the Kalman Filter
T=rows(Y);
N=cols(Y);
ns=cols(beta0);
F=eye(ns);
mu=0;
beta_tt=zeros(T,ns);          %will hold the filtered state variable
ptt=zeros(T,ns,ns);    % will hold its variance
beta11=beta0;
p11=P00;

A=chofac(N,a1');
iA=A\eye(cols(A));

% %%%%%%%%%%%Step 2b run Kalman Filter

for i=1:T
   x=kron(eye(N),X(i,:));


H=diag(hlast(i+1,:));
R=iA*H*iA';
    %Prediction
beta10=mu+beta11*F';
p10=F*p11*F'+Q;
yhat=(x*(beta10)')';                                               
eta=Y(i,:)-yhat;
feta=(x*p10*x')+R;
%updating
ifeta=feta\eye(cols(feta));
K=(p10*x')*ifeta;
beta11=(beta10'+K*eta')';
p11=p10-K*(x*p10);

ptt(i,:,:)=p11;
beta_tt(i,:)=beta11;

end


%%%%%%%%%%%end of Kalman Filter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%step 2c Backward recursion to calculate the mean and variance of the distribution of the state
%vector
chck=-1;
while chck<0
beta2 = zeros(1,ns);   %this will hold the draw of the state variable
wa=randn(1,ns);
roots=zeros(1,1);

i=T;  %period t
p00=squeeze(ptt(i,:,:)); 
beta2(1,:)=beta_tt(i:i,:)+(wa(1:1,:)*cholx(p00));   %draw for beta in period t from N(beta_tt,ptt)
error=Y-X*reshape(beta2(1,:),N*L+2,N);  %var residuals
roots(1)=stability(beta2(1,:)',N,L,2);



if sum(roots)==0
    chck=1;
end
end