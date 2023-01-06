function ldens = lpdfig(x,a,b)
% log INVERSE GAMMA 
% f(x;a,b) = b^a/G(a) x^(-a-1) exp(-b/x)
% where G(.) is the gamma function
% a = mu^2/sigma^2 + 2
% b = mu(a-1)
% philip.liu@anu.edu.au [1/3/2008]

ldens = a*log(b) - log(gamma(a)) + (-a-1).*log(x) - b./x;