
function [bdraw,Q]=getAR(Y,X,B0,Sigma0,sigma2,v0,t0)
V=invpd(invpd(Sigma0)+(1/sigma2)*(X'*X));
M=V*(invpd(Sigma0)*B0+(1/sigma2)*X'*Y); 
	
	chck=-1;
	while chck<0
		bdraw=M+(randn(1,cols(X))*chol(V))';
		if abs(bdraw(1))<=1;
			chck=1;
        end
    end
	res=Y-X*bdraw;
	Q=IG(t0,v0,res);
	