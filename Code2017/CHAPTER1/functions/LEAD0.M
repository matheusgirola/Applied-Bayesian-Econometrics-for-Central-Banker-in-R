function out=lead0(x,p)


[R,C]=size(x);

%Take the first R-p rows of matrix x
x1=x(1:(R-p),:);
x1=x(p+1:end,:);
%Preceed them with p rows of zeros and return
out=[x1;zeros(p,C)];


