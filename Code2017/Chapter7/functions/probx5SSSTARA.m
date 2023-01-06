function accept=probx5SSSTARA(y,x,varcoef1,varcoef2,iamat1,htrial,Smat,e2)
	
   
    
    x1=x;
    x2=x.*repmat(e2,1,cols(x));
    
    
	res2=y-x1*varcoef1-x2*varcoef2;
    
	
A1=iamat1;

sigma1=A1*diag(htrial(1,1).*Smat)*A1';


if rcond(sigma1)<1e-12 
	accept=-1000000;
else
    isigma1=invpd(sigma1);
    
ldet1=logdet(sigma1);

if isnan(ldet1) || isinf(ldet1) || ~isreal(ldet1) ;
 accept=-1000000;
else
w2=(-0.5*ldet1-0.5*(res2*isigma1*res2')); 
accept=(w2);
if isnan(accept) || isinf(accept) || ~isreal(accept);
	accept=-1000000;
end;
end;
end;
