function outx=remSC(x,rho)
L=cols(rho);
out=0;
for j=1:L
    out=out+lag0(x,j).*repmat(rho(:,j),1,cols(x));
end
outx=x-out;