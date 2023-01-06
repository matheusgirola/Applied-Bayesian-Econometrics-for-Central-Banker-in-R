function  f = loglik(b,sigma,y,x)
n=length(b);
z=[];
for i=1:n
e=y-x*b(i);
z=[z;e'*e];
end

f = exp(-0.5*(1/sigma)*z);%*100000000000000*100000000;


