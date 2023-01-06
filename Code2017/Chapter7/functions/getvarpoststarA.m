function [ post,lik1,prior ] = getvarpoststarA( y,x,beta1,beta2,Sbig,iamat,L,param,tarmean,tarvariance,gmean,gvariance,Ystar,ncrit,lamda,EX )
T=rows(y);
N=cols(y);
tar=param(1);
gam=param(2);
LSTAR=1./(1+exp(-gam.*(Ystar-tar)));
e1=(1-LSTAR);
e2=LSTAR;
if sum(e1)<ncrit || sum(e2)<ncrit
    post=-inf;
    lik1=-inf;
    lik2=-inf;
   prior=-inf;
else


lamday=repmat(sqrt(lamda),1,N);
lamdax=repmat(sqrt(lamda),1,cols(x));
Y=y./lamday;
X=x./lamdax;



X1=X;
X2=X.*repmat(e2,1,cols(X));


    sigma1=iamat*diag(Sbig)*iamat';

    lik1=loglikx(reshape(beta1,N*L+EX,N),reshape(beta2,N*L+EX,N),sigma1,Y,X1,X2);
    
    %evaluate prior for the threshold
prior=multivariatenormal(tar,tarmean,tarvariance);
b =gvariance^2/gmean;
a = gmean/b;
priorg= lpdfgam(gam,a,b);
post=lik1+prior+priorg;
if isinf(post) || ~isreal(post) || isnan(post)
    post=-inf;
end
end



