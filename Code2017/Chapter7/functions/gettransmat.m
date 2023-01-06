function [Fx,Mux] = gettransmat( NC,nfact,L,beta2w,beta2e,beta2c,N,EX )

%transition equation
Fx=zeros(((NC+2)*nfact)*L,((NC+2)*nfact)*L);
Mux=zeros(((NC+2)*nfact)*L,1);

%world factor
beta3=reshape(beta2w,N*L+EX,N);
exo3=beta3(N*L+1:end,:);
Mux(1:nfact,:)=exo3';
beta4=beta3(1:N*L,:)';
Fxi=zeros(nfact,((NC+2)*nfact)*L);
jjx=1;
kkx=1;
for jx=1:L
    tmp=beta4(:,jjx:jjx+N-1);
    Fxi(:,kkx:kkx+N-1)=tmp;
    jjx=jjx+N;
    kkx=kkx+(nfact*(NC+2));
end 
Fx(1:nfact,:)=Fxi;
%EA Factor
beta3=reshape(beta2e,N*L+EX,N);
exo3=beta3(N*L+1:end,:);
Mux(nfact+1:nfact*2,:)=exo3';
beta4=beta3(1:N*L,:)';
Fxi=zeros(nfact,((NC+2)*nfact)*L);
jjx=1;
kkx=nfact+1;
for jx=1:L
    tmp=beta4(:,jjx:jjx+N-1);
    Fxi(:,kkx:kkx+N-1)=tmp;
    jjx=jjx+N;
    kkx=kkx+(nfact*(NC+2));
end 
Fx(nfact+1:nfact*2,:)=Fxi;
%country factors

jj=1;
kk=nfact+nfact+1;
for j=1:rows(beta2c)
beta3=reshape(beta2c(j,:),N*L+EX,N);

         
exo3=beta3(N*L+1:end,:);
beta4=beta3(1:N*L,:)';

Fxi=zeros(nfact,((NC+2)*nfact)*L);
jjx=1;
kkx=kk;
for jx=1:L
    tmp=beta4(:,jjx:jjx+N-1);
    Fxi(:,kkx:kkx+N-1)=tmp;
    jjx=jjx+N;
    kkx=kkx+(nfact*(NC+2));
end 


Fx(kk:kk+N-1,:)=Fxi;
Mux(kk:kk+N-1,:)=exo3';
kk=kk+N;
jj=jj+nfact;
end
Fx(((NC+2)*nfact)+1:end,1:end)=eye(((NC+2)*nfact)*(L-1),((NC+2)*nfact)*L);

end

