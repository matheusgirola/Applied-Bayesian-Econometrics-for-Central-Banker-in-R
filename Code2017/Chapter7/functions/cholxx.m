function out=cholxx(x,n)
	temp=x(1:n,1:n);
    
    ctemp=cholx(temp);
	out=zeros(rows(x),cols(x));
	out(1:n,1:n)=ctemp;