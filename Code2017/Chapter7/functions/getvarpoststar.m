function [ post,lik1,lik2,prior ] = getvarpoststar( y,x,beta1,beta2,Sbig1,Sbig2,iamat1,iamat2,L,param,tarmean,tarvariance,gmean,gvariance,Ystar,ncrit,lamda,EX )
T=rows(y);
N=cols(y);
tar=param(1);
gam=param(2);
LSTAR=1./(1+exp(-gam.*(Ystar-tar)));
e1=(1-LSTAR);
e2=LSTAR;



lamday=repmat(sqrt(lamda),1,N);
lamdax=repmat(sqrt(lamda),1,cols(x));
Y=y./lamday;
X=x./lamdax;

Y1=Y.*repmat(e1,1,N);
Y2=Y.*repmat(e2,1,N);

X1=X.*repmat(e1,1,cols(X));
X2=X.*repmat(e2,1,cols(X));


    sigma1=iamat1*diag(Sbig1)*iamat1';
    sigma2=iamat2*diag(Sbig2)*iamat2';

    lik1=loglik(reshape(beta1,N*L+EX,N),sigma1,Y1,X1);
    lik2=loglik(reshape(beta2,N*L+EX,N),sigma2,Y2,X2);
    
    %evaluate prior for the threshold
prior=multivariatenormal(tar,tarmean,tarvariance);
b =gvariance^2/gmean;
a = gmean/b;
priorg= lpdfgam(gam,a,b);
post=lik1+lik2+prior+priorg;
if isinf(post) || ~isreal(post) || isnan(post)
    post=-inf;
end




