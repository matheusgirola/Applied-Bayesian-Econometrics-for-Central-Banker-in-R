function [ post,lik1,lik2,prior ] = getvarpost( Y,X,beta1,beta2,sigma1,sigma2,L,tar,tarmean,tarvariance,Ystar,ncrit )
T=rows(Y);
N=cols(Y);

e1=Ystar<=tar;
e2=Ystar>tar;
if sum(e1) <ncrit ||sum(1-e1)<ncrit
    post=-inf;
    lik=-inf;
    prior=-inf;
else
    
    Y1=Y(e1,:);
    X1=X(e1,:);
    
    Y2=Y(e2,:);
    X2=X(e2,:);
    
   

    lik1=loglik(reshape(beta1,N*L+1,N),sigma1,Y1,X1);
    lik2=loglik(reshape(beta2,N*L+1,N),sigma2,Y2,X2);
    
    %evaluate prior for the threshold
prior=multivariatenormal(tar,tarmean,tarvariance);

post=lik1+lik2+prior;
if isinf(post) || ~isreal(post) || isnan(post)
    post=-inf;
end


end

