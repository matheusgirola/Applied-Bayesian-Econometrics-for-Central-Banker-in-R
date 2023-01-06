function out= drchrnd(a)
m=1;
n=length(a);
r=zeros(n,1);
out=zeros(n,1);
for i=1:n
    r(i)=gamrnd(a(i),m,1,1);
end

out(1:n-1)=r(1:n-1)./sum(r(1:n));
out(n)=1-sum(out(1:n-1));