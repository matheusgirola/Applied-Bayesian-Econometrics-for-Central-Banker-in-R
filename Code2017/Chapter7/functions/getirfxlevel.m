function [ irfw,irfe,irfcc,hlastw,hlaste,hlastc ] =...
    getirfxlevel( fload,varcoef,F,Q,varcoefe,Fe,Qe,varcoefC,FC,QC,Horz,L,Lh,NC ,nfact,sigmaw,sigmac,sigmae,pos)

NN=rows(fload);
T=Horz+L;
ymatw=zeros(T,nfact);
ymate=zeros(T,nfact);
ymatc=zeros(T,nfact);
hlastc=zeros(T,1);
hlaste=zeros(T,1);
hlastw=zeros(T,1);
shockw=zeros(T,nfact);
shockw(L+1,pos)=1;
A0w=chol(sigmaw);
for j=L+1:Horz+L
    hlastw(j,:)=hlastw(j-1,:)*F;
    xhatw=[];
    for i=1:L
        xhatw=[xhatw ymatw(j-i,:)];
    end
    for i=1:Lh
        xhatw=[xhatw hlastw(j-i,:)];
    end
    
    ymatw(j,:)=[xhatw 0]*varcoef+shockw(j,:)*A0w;
end
irfw=[ymatw zeros(T,nfact) zeros(T,nfact*NC)]*fload';
irfw=irfw(L+1:end,:);

%EA factor level shock
shocke=zeros(T,nfact);
shocke(L+1,pos)=1;
A0e=chol(sigmae);
for j=L+1:Horz+L
    hlaste(j,:)=hlaste(j-1,:)*Fe;
    xhate=[];
    for i=1:L
        xhate=[xhate ymate(j-i,:)];
    end
    for i=1:Lh
        xhate=[xhate hlaste(j-i,:)];
    end
    ymate(j,:)=[xhate 0]*varcoefe+shocke(j,:)*A0e;
end
irfe=[ zeros(T,nfact) ymate zeros(T,nfact*NC)]*fload';
irfe=irfe(L+1:end,:);

%country uncertainty shock
irfcc=cell(NC,1);
mm=1;
for m=1:NC
    varcoefc=varcoefC{m};
    Fc=FC{m};
    Qc=QC{m};
    Sc=sigmac{m};
    A0c=chol(Sc);
    
shockc=zeros(T,nfact);
shockc(L+1,pos)=1;
for j=L+1:Horz+L
    hlastc(j,:)=hlastc(j-1,:)*Fc;
    xhatc=[];
    for i=1:L
        xhatc=[xhatc ymatc(j-i,:)];
    end
    for i=1:Lh
        xhatc=[xhatc hlastc(j-i,:)];
    end
    ymatc(j,:)=[xhatc 0]*varcoefc+shockc(j,:)*A0c;
end
   ymatcc=zeros(T,NC*nfact);
   ymatcc(:,mm:mm+nfact-1)=ymatc;
irfc=[ zeros(T,nfact*2) ymatcc]*fload';
irfc=irfc(L+1:end,:);
irfcc{m,1}=irfc;
mm=mm+nfact;
end

