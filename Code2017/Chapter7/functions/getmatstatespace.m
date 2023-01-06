function [ H,Fx,Mux,H0,H1 ] = getmatstatespace( fload,rho,beta2w,beta2c,index,NN,NC,nfact,L,N )
%build matrices of the state-space
id=unique(index);
H0=fload(:,1:nfact);
H1=zeros(NN,NC*nfact);
jj=1;
kk=1;
for j=1:length(id);
    floadc=fload(index==id(j),nfact+1:end);
    H1(kk:kk+rows(floadc)-1,jj:jj+nfact-1)=floadc;
    kk=kk+rows(floadc);
    jj=jj+nfact;
end
H2=[H0 H1];
H3=H2.*-repmat(rho,1,cols(H2));
H=zeros(NN,((NC+1)*nfact)*L);
H(:,1:((NC+1)*nfact)*2)=[H2 H3];

%transition equation
Fx=zeros(((NC+1)*nfact)*L,((NC+1)*nfact)*L);
Mux=zeros(((NC+1)*nfact)*L,1);

%world factor
beta3=reshape(beta2w,N*L+1,N);
mu3=beta3(end,:)';
beta4=beta3(1:N*L,:)';
Fxi=zeros(nfact,((NC+1)*nfact)*L);
jjx=1;
kkx=1;
for jx=1:L
    tmp=beta4(:,jjx:jjx+N-1);
    Fxi(:,kkx:kkx+N-1)=tmp;
    jjx=jjx+N;
    kkx=kkx+(nfact*(NC+1));
end 
Fx(1:nfact,:)=Fxi;
Mux(1:nfact,1)=mu3;


%country factors

jj=1;
kk=nfact+1;
for j=1:rows(beta2c)
beta3=reshape(beta2c(j,:),N*L+1,N);
mu3=beta3(end,:)';
beta4=beta3(1:N*L,:)';

Fxi=zeros(nfact,((NC+1)*nfact)*L);
jjx=1;
kkx=kk;
for jx=1:L
    tmp=beta4(:,jjx:jjx+N-1);
    Fxi(:,kkx:kkx+N-1)=tmp;
    jjx=jjx+N;
    kkx=kkx+(nfact*(NC+1));
end 


Fx(kk:kk+N-1,:)=Fxi;
Mux(kk:kk+N-1,1)=mu3;
kk=kk+N;
end
Fx(((NC+1)*nfact)+1:end,1:end)=eye(((NC+1)*nfact),((NC+1)*nfact)*L);

end

