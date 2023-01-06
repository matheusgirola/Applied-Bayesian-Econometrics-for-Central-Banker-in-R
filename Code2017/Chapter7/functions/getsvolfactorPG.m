function [ hnew,xfix ] = getsvolfactorPG( xfix,fbig,Qbig,mubig,beta2,Sbig,...
    y,x,N,L,EX,mu0,iamat,npart)


NS=cols(xfix);
T=rows(x);

%%%%%%%%%%%%Step 5 Stochastic vol%%%%%%%%
Fmat=zeros(NS,NS);
Fmat(1:1,1:1)=fbig(1);
Qmat=zeros(NS,NS);
Qmat(1:1,1:1)=Qbig;
iQmat=zeros(NS,NS);
iQmat(1:1,1:1)=1/(Qbig);
mumat=zeros(NS,1);
mumat(1,1)=mubig(1,1);
if NS>1
Fmat(2:NS,1:NS-1)=eye(NS-1);
end
cQmat=sqrt(Qmat);
varcoef=reshape(beta2,N*L+EX,N); 

[betax,BB,W]=pftest3F(Fmat,mumat',cQmat,iQmat,iamat,xfix,y,x,npart,log(mu0),varcoef,Sbig);
index=discretesample(W(:,end),1);
xfix=betax(:,index);
hnew=exp(xfix);