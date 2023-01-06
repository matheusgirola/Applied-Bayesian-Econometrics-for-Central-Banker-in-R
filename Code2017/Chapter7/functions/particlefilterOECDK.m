function [ lik,statesw,statese,statesc,statesr,statesp,problem,dens,densexp,likvec,ESS ] = ...
    particlefilterOECDK( y,Qin,fmatw,mumatw,qmatw,Sbigw,...
    fmate,mumate,qmate,Sbige,...
    fmatc,mumatc,qmatc,Sbigc,fmatr,mumatr,qmatr,...
    npart,T,N,L,EX,B00w,P00w,B00e,P00e,B00c,P00c,B00r,P00r,B00p,P00p,p0p,B00f,P00f,nfact,...
    Fx,MUx,iamatw,iamatc,iamate,NC,rho,index,id,indexe)
%matrices of the state space
problem=0;
NSw=cols(P00w);
Fw=zeros(NSw,NSw);
Fw(1:1,1:1)=fmatw(1);
Qw=zeros(NSw,NSw);
Qw(1,1)=qmatw;
cQw=zeros(NSw,NSw);
cQw(1,1)=sqrt(qmatw);
MUw=zeros(NSw,1);
MUw(1)=mumatw;

%EA
NSe=cols(P00e);
Fe=zeros(NSe,NSe);
Fe(1:1,1:1)=fmate(1);
Qe=zeros(NSe,NSe);
Qe(1,1)=qmate;
cQe=zeros(NSe,NSe);
cQe(1,1)=sqrt(qmate);
MUe=zeros(NSe,1);
MUe(1)=mumate;







NSc=cols(P00c);
Fc=diag(fmatc);
Qc=diag(qmatc);

cQc=diag(sqrt(qmatc));
MUc=mumatc;


NN=cols(y);

ns=cols(P00p);
NSf=size(Qin,2);

lik=0;
likvec=zeros(T,1);
ESS=zeros(T,2);

statesw=zeros(T,NSw);
statese=zeros(T,NSe);
statesc=zeros(T,NSc);
statesr=zeros(T,NN);
statesp=zeros(T,ns);
%initial state
b0f=zeros(NN,NSf,npart);
for ii=1:NN
    if indexe(ii)==1
    b0f(ii,:,:)=(repmat(B00f(ii,:),npart,1)+randn(npart,NSf)*P00f)';
    else
    b0f(ii,:,:)=(repmat(B00f(ii,:),npart,1)+[randn(npart,nfact) zeros(npart,nfact) randn(npart,nfact)] *P00f)';
    end
end
b0w=repmat(B00w,npart,1)+[randn(npart,1) ]*(P00w);
b0e=repmat(B00e,npart,1)+[randn(npart,1) ]*(P00e);

b0c=repmat(B00c,npart,1)+[randn(npart,NSc) ]*(P00c);

b0r=repmat(B00r,npart,1)+randn(npart,NN)*P00r;
b0p=repmat(B00p,npart,1)+[randn(npart,nfact*(NC+2)) zeros(npart,ns-(nfact*(NC+2)))]*(P00p);

for i=1:T
 
    dens=zeros(npart,1);
    partr=randn(npart,NN);
    parthw=zeros(npart,NSw);
   
    parthw(:,1)=randn(npart,1);
    parthw=parthw*cQw;
    
    parthe(:,1)=randn(npart,1);
    parthe=parthe*cQe;
    
    parthc=randn(npart,NSc);
    parthc=parthc*cQc;
    
    parthf=zeros(NN,NSf,npart);
    for ii=1:NN
        Qi=squeeze(Qin(ii,:,:));
        if indexe(ii)==1
            cQi=cholx(Qi);
            tmp=randn(npart,NSf)*cQi;
        else
            cQi=cholx(Qi(1:nfact*2,1:nfact*2));
            tmp1=randn(npart,NSf-nfact)*cQi;
            tmp=[tmp1(:,1:nfact) zeros(npart,nfact) tmp1(:,nfact+1:end)];
        end
        parthf(ii,:,:)=tmp';
    end

    
    
parfor j=1:npart
 %run the prediction step of the Kalman filter    
    dens(j)=predictionstep(y(i,:),b0f(:,:,j),Fx,MUx,b0p(j,:),squeeze(p0p(j,:,:)),...
    b0r(j,:),b0c(j,:), b0w(j,:),b0e(j,:),iamatc,iamatw,iamate,Sbigc,Sbigw,Sbige,...
    NC,nfact,rho,index,id,L);





end
%calculate weights and re-sample
  [densexp,sfactor]=safeexp(dens);
    sdens=sum(densexp);
    if sdens==0 || isinf(sdens)
        lik=-999;
        problem=1;
        return
    end
        
   prob=densexp./sdens;
   
   index1=discretesample(prob,npart);
    tmp1=length(unique(index1));
    tmp2=1./sum(prob.^2);
    ESS(i,1)=tmp2;
    ESS(i,2)=tmp1;
  
   %re-sample
   b1f=b0f(:,:,index1);
   b1w=b0w(index1,:);
   b1e=b0e(index1,:);
    b1c=b0c(index1,:);
   b1r=b0r(index1,:);
   b1p=b0p(index1,:);
   p1p=p0p(index1,:,:);
   
   parfor j=1:npart
       
       
       [b0w(j,:),b0c(j,:),b0r(j,:),b0p(j,:),p0p(j,:,:),b0f(:,:,j)]=...
           getdensOECDk(b1w(j,:),parthw(j,:),b1e(j,:),parthe(j,:),...
    b1c(j,:),parthc(j,:),b1r(j,:),partr(j,:),...
    Sbigw,Sbige,Sbigc,b1p(j,:),squeeze(p1p(j,:,:)),y(i,:),MUw,Fw,MUe,Fe,MUc,...
    Fc,mumatr,fmatr,qmatr,...
    nfact,Fx,MUx,iamatw,iamate,iamatc,NC,b1f(:,:,j),parthf(:,:,j),NN,NSf,rho,index,id,L);

% [hneww,hnewc,rnew,pmatnew,p0new,floadnew]= getdensOECDk(b0w,parthw,b0e,parthe,...
%     b0c,parthc,b0r,partr,...
%     Sbigw,Sbige,Sbigc,b0p,p0p,y,MUw,Fw,MUe,Fe,MUc,Fc,MUr,Fr,Qr,...
%     nfact,Fx,Mux,iAAw,iAAe,iAAc,NC,b0f,parthf,NN,NSf,rho,index,id,L)

   end
      

  
   
   
   liki=sdens/npart;
      num2str([i (log(liki)+sfactor) sfactor])

lik=lik+(log(liki)+sfactor);
likvec(i)=(log(liki)+sfactor);
statesw(i,:)=sum(b0w.*repmat(prob,1,NSw))';
statese(i,:)=sum(b0e.*repmat(prob,1,NSe))';

statesc(i,:)=sum(b0c.*repmat(prob,1,NSc))';
statesr(i,:)=sum(b0r.*repmat(prob,1,cols(statesr)))';
statesp(i,:)=sum(b0p.*repmat(prob,1,cols(statesp)))';
end

