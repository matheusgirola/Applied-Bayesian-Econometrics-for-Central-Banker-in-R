function accept=probx5SSM(y,x,varcoef,htrial)
	
	

 w2=-0.5*log(htrial(1)) - ((y-x*varcoef)^2)/(2*htrial(1)); 
accept=(w2);
if isnan(accept) || isinf(accept) || ~isreal(accept);
	accept=-1000000;
end;

