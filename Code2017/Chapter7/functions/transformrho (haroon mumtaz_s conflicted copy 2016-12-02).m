function out=transformrho(x,rho)
L=cols(rho);
rhs=0;
for i=1:L
    xl=lag0(x,i).*rho(i);
    rhs=rhs+xl;
end
% size(rhs)
out=x-rhs;