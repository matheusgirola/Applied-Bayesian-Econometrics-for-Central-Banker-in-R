
function [bdraw,Q]=getARmc(Y,X,B0,Sigma0,sigma2,v0,t0)
V=invpd(invpd(Sigma0)+(1/sigma2)*(X'*X));
M=V*(invpd(Sigma0)*B0+(1/sigma2)*X'*Y); 
	
	reps=100;	
	chck=-1;
    trys=1;
	while chck<0 &&trys <reps
		bdraw=M+(randn(1,cols(X))*chol(V))';
        S=abs(bdraw(1))>=1;

		if S==0;
			chck=1;
        else
            trys=trys+1;
        end
    end
    
    if chck<0;
       bdraw=B0;
    end
	res=Y-X*bdraw;
	Q=IG(t0,v0,res);
	