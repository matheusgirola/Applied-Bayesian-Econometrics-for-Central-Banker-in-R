function  H = getH( fload,nfact,NN,NC,index,id,rho,L )
%build matrices of the state-space (hard coded for serial correlation order
%1 in idiosyncratic errors)
H0=fload(:,1:nfact);
H0e=fload(:,nfact+1:nfact*2);
H1=zeros(NN,NC*nfact);
jj=1;
kk=1;
for j=1:length(id);
    floadc=fload(index==id(j),(nfact*2)+1:end);
    H1(kk:kk+rows(floadc)-1,jj:jj+nfact-1)=floadc;
    kk=kk+rows(floadc);
    jj=jj+nfact;
end
H2=[H0 H0e H1];
H3=H2.*-repmat(rho,1,cols(H2)); %this would change for longer lags in serial correlation
H=zeros(NN,((NC+2)*nfact)*L);
H(:,1:((NC+2)*nfact)*2)=[H2 H3];

end

