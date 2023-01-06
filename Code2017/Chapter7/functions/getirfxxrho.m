function [ irf ] = getirfxx( fload,F,Q,Horz,L,Lh )
NN=rows(fload);
nfact=cols(fload);
yhat=zeros(Horz+L,nfact);
hhat=zeros(Horz+L,1);
shock=zeros(Horz+L,1);
shock(L+1)=1;
for j=L+1:Horz+L
    hhat(j,:)=hhat(j-1,:)*F+shock(j,:)*sqrt(Q);
    
end
irf=hhat.*fload;

irf=irf(L+1:end,:);

end

