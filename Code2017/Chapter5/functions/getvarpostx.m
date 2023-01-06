function [ post,lik1,prior ] = getvarpostx( Y,X,beta1,beta2,sigma1,L,param,tarmean,tarvariance,gmean,gvariance,Ystar,ncrit )
T=rows(Y);
N=cols(Y);

tar=param(1);
gam=param(2);

LSTAR=1./(1+exp(-gam.*(Ystar-tar)));
    e1=(1-LSTAR);
    e2=LSTAR;
    
    
X1=X.*repmat(e1,1,cols(X));
    
   
    X2=X.*repmat(e2,1,cols(X));
    
   

    lik1=loglikx(reshape(beta1,N*L+1,N),reshape(beta2,N*L+1,N),sigma1,Y,X1,X2);
    
    
    %evaluate prior for the threshold
prior=multivariatenormal(tar,tarmean,tarvariance);
b =gvariance^2/gmean;
a = gmean/b;
priorg= lpdfgam(gam,a,b);

post=lik1+prior+priorg;
if isinf(post) || ~isreal(post) || isnan(post)
    post=-inf;
end




