function [ irf ] = getirfxxL( fload,F,Q,Horz,L,Lh )
NN=rows(fload);
nfact=cols(fload);
yhat=zeros(Horz+L,1);
hhat=zeros(Horz+L,1);
shock=zeros(Horz+L,1);
shock(L+1)=1;
for j=L+1:Horz+L
    hhat(j,:)=hhat(j-1,:)*F+shock(j,:)*sqrt(Q);
    xhat=[];
    for i=0:Lh
        xhat=[xhat hhat(j-i,:)];
    end
    yhat(j,:)=xhat*fload';
end


irf=yhat(L+1:end,:);

end

