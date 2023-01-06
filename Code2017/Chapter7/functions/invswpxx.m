function out=invswpxx(x,n)
	temp=x(1:n,1:n);
   itemp=invpd(temp);
	out=zeros(rows(x),cols(x));
	out(1:n,1:n)=itemp;