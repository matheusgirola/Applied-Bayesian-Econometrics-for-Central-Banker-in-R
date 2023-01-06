
function bdraw=getARxx(Y,X,B0,Sigma0,sigma2)
problem=0;
V=invpd(invpd(Sigma0)+(1/sigma2)*(X'*X));
M=V*(invpd(Sigma0)*B0+(1/sigma2)*X'*Y); 

	chck=-1;
    trys=1;
	while chck<0
		bdraw=M+(randn(1,cols(X))*chol(V))';
%         S=stability(bdraw,1,cols(X),0);
S=abs(bdraw(1))>=1;
		if S==0;
			chck=1;
        end
    end
    
    
	