function [yd,xd,b0,s0]=...
    getdums(LAMDAPEVIEWS,...
    TAUPEVIEWS,L)
N=1;

lamdaP=LAMDAPEVIEWS;  %This controls the tightness of the priors on the first lag
tauP=TAUPEVIEWS;  % this controls the tightness of the priors on sum of coefficients
muP=0;
sigmaP=ones(N,1);
deltaP=ones(N,1).*0.8;

%dummy data to implement priors see http://ideas.repec.org/p/ecb/ecbwps/20080966.html
[yd,xd] = create_dum(lamdaP,tauP,deltaP,L,muP,sigmaP,N);
% size(yd)
% size(xd)
b0=vec(xd\yd);

  b01=reshape(b0,N*L,N);
  e0=yd-xd*b01;
%   S=e0'*e0;
S=diag(sigmaP.^2);
  s0=kron(S,pinv(xd'*xd));