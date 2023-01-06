function out = multivariatenormal(x,mu,sigma)
k=rows(x);
res=x-mu;
isigma=invpd(sigma);
expterm=-0.5*res'*isigma*res;
constant=log(1/(2*(pi^k/2)))-0.5*logdet(sigma);
out=constant+expterm;


