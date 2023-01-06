function  f = postlik(b,sigma,y,x,X,m,s)
n=length(b);
z=[];
for i=1:n
e=y-x*b(i);
z=[z;e'*e];
end
f1=-((X-m)./s).^2/2;
%-log(s)-log(2*pi)/2
f = exp(-0.5*(1/sigma)*z);



f=f+exp(f1)';

