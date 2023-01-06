function accept=probx5SSSTAR(y,x,varcoef1,varcoef2,iamat1,iamat2,htrial,Smat,e2)
	
    e1=1-e2;
    y1=y.*repmat(e1,1,cols(y));
    y2=y.*repmat(e2,1,cols(y));
    x1=x.*repmat(e1,1,cols(x));
    x2=x.*repmat(e2,1,cols(x));
    
    
	res21=y1-x1*varcoef1;
    res22=y2-x2*varcoef2;
	
A1=iamat1;
A2=iamat2;

sigma1=A1*diag(htrial(1,1).*Smat)*A1';
sigma2=A2*diag(htrial(1,1).*Smat)*A2';


if rcond(sigma1)<1e-12 || rcond(sigma2)<1e-12
	accept=-1000000;
else
    isigma1=invpd(sigma1);
    isigma2=invpd(sigma2);
ldet1=logdet(sigma1);
ldet2=logdet(sigma2);
if isnan(ldet1) || isinf(ldet1) || ~isreal(ldet1) || isnan(ldet2) || isinf(ldet2) || ~isreal(ldet2);
 accept=-1000000;
else
w2=(-0.5*ldet1-0.5*(res21*isigma1*res21'))+...
    (-0.5*ldet2-0.5*(res22*isigma2*res22'));  %numerator
 
accept=(w2);
if isnan(accept) || isinf(accept) || ~isreal(accept);
	accept=-1000000;
end;
end;
end;
