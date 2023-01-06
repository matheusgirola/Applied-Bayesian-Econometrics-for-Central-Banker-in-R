function[beta_tt,ptt]=kalmafiltervar(Y,X,Q,amat,hlast,beta0,P00)
T=rows(Y);
N=cols(Y);
%%Step 2a Set up matrices for the Kalman Filter

ns=cols(beta0);
F=eye(ns);
mu=0;
beta_tt=zeros(T,ns);          %will hold the filtered state variable
ptt=zeros(T,ns,ns);    % will hold its variance
beta11=beta0;
p11=P00;


% %%%%%%%%%%%Step 2b run Kalman Filter

for i=1:T
   x=kron(eye(N),X(i,:));

a=amat(i,:);
A=chofac(N,a');
H=diag(hlast(i+1,:));
R=invpd(A)*H*invpd(A)';
    %Prediction
beta10=mu+beta11*F';
p10=F*p11*F'+Q;
yhat=(x*(beta10)')';                                               
eta=Y(i,:)-yhat;
feta=(x*p10*x')+R;
%updating
K=(p10*x')*invpd(feta);
beta11=(beta10'+K*eta')';
p11=p10-K*(x*p10);

ptt(i,:,:)=p11;
beta_tt(i,:)=beta11;

end

