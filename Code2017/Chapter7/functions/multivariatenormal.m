function out = multivariatenormal(x,mu,sigma)
k=rows(x);
res=x-mu;
isigma=pinv(sigma);
expterm=-0.5*res'*isigma*res;
constant=- k*log(2*pi)/2-0.5*logdet(sigma);%log(det(sigma));
out=constant+expterm;


