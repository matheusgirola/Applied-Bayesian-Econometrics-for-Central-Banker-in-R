function [ lik,statesw,statesc,statesr,statesp,problem,dens,densexp,likvec,ESS ] = ...
    particlefilterOECDK( y,H,fmatw,mumatw,qmatw,Sbigw,...
    fmatc,mumatc,qmatc,Sbigc,fmatr,mumatr,qmatr,...
    npart,T,N,L,EX,B00w,P00w,B00c,P00c,B00r,P00r,B00p,P00p,p0p,nfact,...
    Fx,MUx,iamatw,iamatc,NC)
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

NSc=cols(P00c);
Fc=diag(fmatc);
Qc=diag(qmatc);

cQc=diag(sqrt(qmatc));
MUc=mumatc;


NN=cols(y);

ns=cols(P00p);

lik=0;
likvec=zeros(T,1);
ESS=zeros(T,2);
statesw=zeros(T,NSw);
statesc=zeros(T,NSc);
statesr=zeros(T,NN);
statesp=zeros(T,ns);
%initial state
b0w=repmat(B00w,npart,1)+[randn(npart,1) ]*(P00w);
b0c=repmat(B00c,npart,1)+[randn(npart,NSc) ]*(P00c);

b0r=repmat(B00r,npart,1)+randn(npart,NN)*P00r;
b0p=repmat(B00p,npart,1)+[randn(npart,nfact*(NC+1)) zeros(npart,ns-(nfact*(NC+1)))]*(P00p);

for i=1:T
  
    dens=zeros(npart,1);
%     hneww=zeros(npart,NSw);
%     hnewc=zeros(npart,NSc);
%     rnew=zeros(npart,NN);
%     pmatnew=zeros(npart,ns);
    partr=randn(npart,NN);
    parthw=zeros(npart,NSw);
   
    parthw(:,1)=randn(npart,1);
    parthw=parthw*cQw;
    parthc=randn(npart,NSc);
    parthc=parthc*cQc;
%     partp=zeros(npart,ns);
%     partp(:,1:nfact*(NC+1))=randn(npart,nfact*(NC+1));
parfor j=1:npart
 %run the prediction step of the Kalman filter    
    dens(j)=predictionstep(y(i,:),H,Fx,MUx,b0p(j,:),squeeze(p0p(j,:,:)),...
    b0r(j,:),b0c(j,:), b0w(j,:),iamatc,iamatw,Sbigc,Sbigw,...
    NC,nfact);
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
   
   b1w=b0w(index1,:);
    b1c=b0c(index1,:);
   b1r=b0r(index1,:);
   b1p=b0p(index1,:);
   p1p=p0p(index1,:,:);
   
   parfor j=1:npart
       
       
       [b0w(j,:),b0c(j,:),b0r(j,:),b0p(j,:),p0p(j,:,:)]=getdensOECDk(b1w(j,:),parthw(j,:),...
    b1c(j,:),parthc(j,:),b1r(j,:),partr(j,:),...
    Sbigw,Sbigc,b1p(j,:),squeeze(p1p(j,:,:)),y(i,:),MUw,Fw,MUc,Fc,H,mumatr,fmatr,qmatr,...
    nfact,Fx,MUx,iamatw,iamatc,NC);



   end
      

  
   
   
   
   liki=sdens/npart;
lik=lik+(log(liki)+sfactor);
likvec(i)=(log(liki)+sfactor);
statesw(i,:)=sum(b0w.*repmat(prob,1,NSw))';
statesc(i,:)=sum(b0c.*repmat(prob,1,NSc))';
statesr(i,:)=sum(b0r.*repmat(prob,1,cols(statesr)))';
statesp(i,:)=sum(b0p.*repmat(prob,1,cols(statesp)))';
end

