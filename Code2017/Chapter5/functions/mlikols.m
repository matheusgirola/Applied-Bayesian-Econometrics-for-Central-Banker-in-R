function mlik=mlikols(B0,SIGMA0,V0,D0,Y,X)
%mean of the gamma prior
S0=V0/D0;
%ols quantities
bhat=X\Y;
T=rows(X);
K=cols(X);
v=T-K;
vbar=v+V0;
ehat=Y-X*bhat;
s=(ehat'*ehat)/(v);
%Mlik calculation
vs=(V0*S0)+(v*s)+(bhat-B0)'*inv(SIGMA0+inv(X'*X))*(bhat-B0);
c=(gamma(vbar/2)*(V0*S0)^(V0/2))/(gamma(V0/2)*(pi^(T/2)));
VV=inv(SIGMA0+inv(X'*X));
mlik=c*((det(VV)/det(SIGMA0))^0.5)*((vs)^(-vbar/2));