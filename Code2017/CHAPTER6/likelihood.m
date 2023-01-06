function out=likelihood(theta,y)

N=length(theta);
N1=N-3;  %last 3 elements of theta are standard deviation of shocks
t=size(y,1);

%model solution to obtain state space form
[  PP, QQ, PROBLEM ] = model_solve( theta(1:N1));

Sigma=diag(theta(N1+1:end));

if PROBLEM
    out=-inf;
else


%COMPUTE MATRICES OF THE STATE SPACE
%Y=H*S
%  S[t]=F*S[t-1]+eta
%  Var(eta)=Q
%Q=QQ*Sigma*QQ'

H=zeros(3,6);
H(1,1)=1;
H(2,2)=1;
H(3,3)=1;

F=PP;

Q=QQ*Sigma*QQ';
mu=0;


lik=0;
%filter
beta11=zeros(1,6);
p11=eye(6);
ptt=zeros(t,6,6);

R=0;


for i=1:t
    %Prediction
beta10=mu+beta11*F';
p10=F*p11*F'+Q;
yhat=(H*(beta10)')';                                                
eta=y(i,:)-yhat;
feta=(H*p10*H')+R;
%updating
rc=rcond(feta);
if rc<1e-15
    out=-inf;
    return;
else
    ifeta=inv(feta);
end
K=(p10*H')*ifeta;
beta11=(beta10'+K*eta')';
p11=p10-K*(H*p10);
ptt(i,:,:)=p11;

liki=-0.5*log(2*pi)-0.5*log(det(feta))+(-0.5*(eta)*ifeta*(eta'));

if isreal(liki) && (~isinf(liki))
    lik=lik+liki;
else
    lik=lik-10;
end

end
out=lik;
end

