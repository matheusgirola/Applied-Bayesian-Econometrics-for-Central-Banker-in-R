function accept=probx5SS(y,x,varcoef,iamat,htrial,Smat)
	
	res2=y-x*varcoef;
	
A=iamat;
N=cols(A);
sigma=A*diag(htrial(1,1).*Smat)*A';


if rcond(sigma)<1e-12
	accept=-1000000;
else
    isigma=invpd(sigma);
ldet=logdet(sigma);
if isnan(ldet) || isinf(ldet) || ~isreal(ldet);
 accept=-1000000;
else
w2=-0.5*ldet-0.5*(res2*isigma*res2');  %numerator
 
accept=(w2);
if isnan(accept) || isinf(accept) || ~isreal(accept);
	accept=-1000000;
end;
end;
end;
