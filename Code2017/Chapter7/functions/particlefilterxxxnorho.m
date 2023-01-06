function [ lik,states,statesr,statesp,problem,dens,densexp,likvec ] = ...
    particlefilterxxxnorho( y,H,bload,Sbig,fmat,mumat,qmat,...
    fmatr,mumatr,qmatr,...
    npart,T,N,L,EX,B00,P00,B00r,P00r,B00p,P00p,nfact,...
    varcoef,iamat)

NN=cols(y);
nsr=size(B00r,2);
%matrices of the state space
problem=0;
NS=cols(P00);
F=zeros(NS,NS);
F(1:1,1:1)=fmat(1);
F(2:NS,1:NS-1)=eye(NS-1);

Q=zeros(NS,NS);
Q(1,1)=qmat;
cQ=zeros(NS,NS);

cQ(1,1)=sqrt(qmat);
MU=zeros(NS,1);
MU(1)=mumat;




ns=cols(P00p);

lik=0;
likvec=zeros(T,1);
states=zeros(T,NS);
statesr=zeros(T,NN);
statesp=zeros(T,ns);
%initial state
b0=repmat(B00,npart,1)+[randn(npart,1) zeros(npart,NS-1)]*(P00);
b0r=zeros(npart,NN,nsr);
for ii=1:nsr
    
b0r(:,:,ii)=repmat(B00r(:,ii)',npart,1)+randn(npart,NN)*P00r(ii,ii);
end
b0p=repmat(B00p,npart,1)+[randn(npart,nfact) zeros(npart,ns-nfact)]*(P00p);
for i=1:T
  
    dens=zeros(npart,1);
    hnew=zeros(npart,NS);
    rnew=zeros(npart,NN,nsr);
    pmatnew=zeros(npart,ns);
    partr=randn(npart,NN);
    parth=zeros(npart,NS);
    parth(:,1)=randn(npart,1);
    parth=parth*cQ;
    partp=zeros(npart,ns);
    partp(:,1:N)=randn(npart,nfact);
parfor j=1:npart
     
    

[dens(j),hnew(j,:),rnew(j,:,:),pmatnew(j,:)]= ...
    getdensnorho(b0(j,:),parth(j,:),...
    b0r(j,:,:),partr(j,:),...
    Sbig,b0p(j,:),partp(j,:),N,L,EX,ns,y(i,:),...
    MU,F,H,mumatr,fmatr,qmatr,nfact,varcoef,iamat,bload,nsr,NN);
      
end
    [densexp,sfactor]=safeexp(dens);
    sdens=sum(densexp);
    if sdens==0 || isinf(sdens)
        lik=-999;
        problem=1;
        return
    end
        
   prob=densexp./sdens;
   
   index1=discretesample(prob,npart);
   
   
   %re-sample
   
   b0=hnew(index1,:);
   b0r=rnew(index1,:,:);
   b0p=pmatnew(index1,:);
   
   
   
   liki=sdens/npart;
lik=lik+(log(liki)+sfactor);
likvec(i)=(log(liki)+sfactor);
states(i,:)=sum(b0.*repmat(prob,1,NS))';
statesr(i,:)=sum(squeeze(b0r(:,:,1)).*repmat(prob,1,cols(statesr)))';
statesp(i,:)=sum(b0p.*repmat(prob,1,cols(statesp)))';
end

