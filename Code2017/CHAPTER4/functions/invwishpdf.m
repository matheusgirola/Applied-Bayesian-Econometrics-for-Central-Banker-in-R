function logpdf = invwishpdf(SIGMA,A,V)
% SIGMA = argument matrix, should be positive definite
% A = p x p symmetric, postitive definite "scale" matrix 
% V = "precision" parameter = "degrees of freeedom"


k=size(SIGMA,1);

logexpterm = -.5*trace(inv(SIGMA)*A) ;
logdetAterm = log(det(A))*(V/2) ;
logdetSterm = log(det(SIGMA))*(-(V+k+1)/2) ;
logtwoterm = log(2)*((V*k)/2) ;
logpiterm = log(pi)*(((k-1)*k)/4) ;

sumgamln=sumgamlnx(k,V);

logpdf = logexpterm + logdetSterm +logdetAterm-( logtwoterm + logpiterm + sumgamln  ) ;

