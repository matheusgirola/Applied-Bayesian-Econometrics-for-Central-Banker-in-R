function [ post,lik1,lik2,prior ] = getvarpost( y,x,beta1,beta2,Sbig1,Sbig2,iamat1,iamat2,L,tar,tarmean,tarvariance,Ystar,ncrit,lamda,EX )
T=rows(y);
N=cols(y);
lamday=repmat(sqrt(lamda),1,N);
lamdax=repmat(sqrt(lamda),1,cols(x));
Y=y./lamday;
X=x./lamdax;
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
    sigma1=iamat1*diag(Sbig1)*iamat1';
    sigma2=iamat2*diag(Sbig2)*iamat2';

    lik1=loglik(reshape(beta1,N*L+EX,N),sigma1,Y1,X1);
    lik2=loglik(reshape(beta2,N*L+EX,N),sigma2,Y2,X2);
    
    %evaluate prior for the threshold
prior=multivariatenormal(tar,tarmean,tarvariance);

post=lik1+lik2+prior;
if isinf(post) || ~isreal(post) || isnan(post)
    post=-inf;
end


end

