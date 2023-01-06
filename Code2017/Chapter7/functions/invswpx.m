function out=invswpx(x,n)
	temp=x(1:n,1:n);
    
%  	itemp=pinv(temp)
    itemp=diag(1./diag(temp));
	out=zeros(rows(x),cols(x));
	out(1:n,1:n)=itemp;