function ldens = lpdfig1(x,a,b)
% log INVERSE GAMMA (type 1)
%
% X ~ IG1(s,nu)
% X = sqrt(Y) where Y ~ IG2(s,nu) and Y = inv(Z) with Z ~ G(nu/2,2/s) (Gamma distribution) 
% 
% See L. Bauwens, M. Lubrano and J-F. Richard [1999, appendix A] for more
% details.
%
% stephane.adjemian@cepremap.cnrs.fr [01/16/2004]

% ldens = log(2) - gammaln(nu/2) - (nu/2).*log(2/s) - (nu+1)*log(x) - .5*s./(x.^2);

  ldens = log(2) - gammaln(b/2) + (b/2)*log(b*a^2/2) - ( (b+1)/2 )*log(x^2) - b*a^2/(2*x^2);