clear
addpath('../functions');
dfolder='./data/';
sfolder='./results/';
clc
maxfile=372; 
tic;
%inputs
REPS=55000;%total reps
BURN=50000;%burn in
mreps=1;
Horz=40;
NC=4;
NN=NC*20;  %20 series per country
index=vec(repmat((1:NC),20,1));
id=unique(index);
Qw=eye(3).*0.01;
gw=0.01;
Qc=eye(3).*0.01;
gc=0.01;
Qe=eye(1).*0.01;
ge=0.01;

T0=100;
Tx=220;
nfact=1;
T=Tx+T0;
hlastc=zeros(T,NC);
hlaste=zeros(T,NN);
hlastw=zeros(T,1);
for i=2:T
for mm=1:NC
hlastc(i,mm)=hlastc(i-1,mm)+randn(1,1).*sqrt(gc);
end
hlastw(i)=hlastw(i-1)+randn(1,1).*sqrt(gw);
for mm=1:NN
hlaste(i,mm)=hlaste(i-1,mm)+randn(1,1).*sqrt(ge);
end
end
%TVP coefficients
rhoc=zeros(T,3,NC);
for j=1:NC
    rhoc(1,:,j)=[1.2 -0.3 0];
end
rhow=zeros(T,3);
rhow(1,:)=[1.2 -0.3 0];
rhoe=zeros(T,NN);
for i=2:T
for mm=1:NC
    chck=-1;
    while chck<0
    tmp=rhoc(i-1,:,mm)+randn(1,3)*chol(Qc);
    ee=stability(tmp,1,2,1);
    if ee==0
        chck=1;
    end
    end
    rhoc(i,:,mm)=tmp;
end
end

for i=2:T
for mm=1:NN
    chck=-1;
    while chck<0
    tmp=rhoe(i-1,mm)+randn(1,1)*chol(Qe);
    ee=stability(tmp,1,1,0);
    if ee==0
        chck=1;
    end
    end
    rhoe(i,mm)=tmp;
end
end

for i=2:T

    chck=-1;
    while chck<0
    tmp=rhow(i-1,:)+randn(1,3)*chol(Qw);
    ee=stability(tmp,1,2,1);
    if ee==0
        chck=1;
    end
    end
    rhow(i,:)=tmp;
end
   
    


floadw=randn(NN,nfact).*sqrt(1);
floadc=randn(NN,nfact).*sqrt(0.5);
tmp=[floadw floadc];
tmp1=[];
for k=1:length(id)
  idx=index==id(k);
  tmp2=tmp(idx,:);
  tmp2(1:nfact*2,:)=eye(nfact*2);
  tmp1=[tmp1;tmp2];
end
floadw=tmp1(:,1:nfact);
floadc=tmp1(:,nfact+1:end);


for m=1:mreps


ymatw=zeros(T,nfact);
ymate=zeros(T,NN);
ymatc=zeros(T,nfact*NC);

xmat=zeros(T,NN);

for i=3:T

        
fload=zeros(NN,(NC*nfact)+nfact);
fload(:,1:nfact)=floadw;

kk=(nfact)+1;
ii=1;
for k=1:length(id)
    idx=index==id(k);
    floadx=floadc(idx,:);
    fload(ii:ii+rows(floadx)-1,kk:kk+nfact-1)=floadx;
    ii=ii+rows(floadx);
    kk=kk+nfact;
end
    


sigmaw=exp(hlastw(i));
csigmaw=chol(sigmaw);


mmm=1;
for mm=1:NC
sigmac=exp(hlastc(i,mm));
csigmac=chol(sigmac);
Bc=rhoc(i,:,mm);
ymatc(i,mmm:mmm+nfact-1)=[ymatc(i-1,mmm:mmm+nfact-1) ymatc(i-2,mmm:mmm+nfact-1)   1]*Bc'+randn(1,nfact)*csigmac;
mmm=mmm+nfact;
end
Bw=rhow(i,:);
ymatw(i,:)=[ymatw(i-1,:) ymatw(i-2,:)   1]*Bw'+randn(1,nfact)*csigmaw;
for mm=1:NN
    Be=rhoe(i,mm);
    sigmae=exp(hlaste(i,mm));
    csigmae=chol(sigmae);
ymate(i,mm)=[ymate(i-1,mm) ]*Be'+randn(1,1)*csigmae;
end

xmat(i,:)=[ymatw(i,:)  ymatc(i,:)]*fload'+ymate(i,:);

%variance decomposition
%unconditional variance world factor
volf=volatility( Bw',sigmaw,1,2 );
volfw(i,:)=diag(volf)';

%unconditional variance country factor
mmm=1;
for mm=1:NC
    sigmac=exp(hlastc(i,mm));
    Bc=rhoc(i,:,mm)';
    volcm=volatility( Bc,sigmac,1,2 );
    volct(i,mmm:mmm+nfact-1)=diag(volcm)';
    mmm=mmm+nfact;
end

%unconditional volatility idiosyncratic factor
for mm=1:NN
    sigmae=exp(hlaste(i,mm));
    Be=rhoe(i,mm)';
    volem=volatility( [Be;0],sigmae,1,1 );
    volet(i,mm)=diag(volem)';
end


floadsquared=fload.^2;
totalvol(i,:)=[volfw(i,:) volct(i,:)]*floadsquared'+volet(i,:);
totalvolw(i,:)=volfw(i,:) *floadsquared(:,1:nfact)';


end


ymat=[ymatw ymatc ymate];
hlast=[hlastw hlastc hlaste];
pmat00x=ymat(T0+1:end,:);
dataS=xmat(T0+1:end,:);
hlast0x=hlast(T0+1:end,:);
fload00x=fload;
vdecompwx=totalvolw(T0+1:end,:)./totalvol(T0+1:end,:);

% elast00=elast(T0+1:end,:);
% 
save(strcat(dfolder,'dataxx0',num2str(m)),'dataS','hlast0x','index','pmat00x','vdecompwx');
end


