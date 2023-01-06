function out=transformrho(x,rho)
L=cols(rho);
rhs=0;
for i=1:L
    xl=lag0(x,i).*rho(i);
    rhs=rhs+xl;
end
% size(rhs)
out1=x(L+1:end,:)-rhs(L+1:end,:);

out=[repmat(out1(1,:),L,1);out1];