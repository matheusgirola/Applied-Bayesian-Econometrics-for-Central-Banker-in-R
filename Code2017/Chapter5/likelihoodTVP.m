function [out,beta_tt]=likelihoodTVP(theta,y,x)

%extract parameters of the state space
out=1000000000;
if sum(theta(5:end)<0)==0 && sum(abs(theta(1:2))>1)==0 
F=zeros(2,2);
F(1,1)=theta(1);
F(2,2)=theta(2);

mu=zeros(1,2);
mu(1,1)=theta(3);
mu(1,2)=theta(4); 
    
    
r=(theta(5));
Q=zeros(2,2);
Q(1,1)=(theta(6));
Q(2,2)=(theta(7));

t=rows(y);
lik=0;
%filter
beta0=zeros(1,2);
p00=eye(2);
beta_tt=[];
ptt=zeros(t,2,2);

beta11=beta0;
p11=p00;


for i=1:t
    H=x(i,:);
    %Prediction
beta10=mu+beta11*F';
p10=F*p11*F'+Q;
yhat=(H*(beta10)')';                                                
eta=y(i,:)-yhat;
feta=(H*p10*H')+r;
%updating
K=(p10*H')*inv(feta);
beta11=(beta10'+K*eta')';
p11=p10-K*(H*p10);
beta_tt=[beta_tt;beta11];
ptt(i,:,:)=p11;

liki=-0.5*log(2*pi)-0.5*log(det(feta))+(-0.5*(eta)*inv(feta)*(eta'));

if isreal(liki) & (1-isinf(liki))
    lik=lik+liki;
else
    lik=lik-10;
end

end
out=-lik;
if isnan(out)|| 1-isreal(out) || isinf(out)
    out=1000000000;
end
end

