function bdraw=getregx(Y,X,B0,SIGMA0,sigma2)
isigma0=invpd(SIGMA0);
	isigma2=1/sigma2;
	xx=X'*X;
	V=invpd(isigma0+isigma2*xx);
	M=V*(isigma0*B0+isigma2*X'*Y);
    bdraw=M+(randn(1,cols(X))*chol(V))';
    