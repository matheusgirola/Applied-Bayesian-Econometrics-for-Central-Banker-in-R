function [ lik,states ] = ...
    particlefilter( y,x,varcoef1,varcoef2,iamat1,iamat2,Sbig,fmat,mumat,qmat,...
    npart,T,N,L,EX,B00,P00,e1)
%matrices of the state space
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


H=eye(N,NS);
lik=0;
states=zeros(T,NS);
%initial state
b0=repmat(B00,npart,1)+[randn(npart,1) zeros(npart,NS-1)]*chol(P00);
for i=1:T
    
    dens=zeros(npart,1);
	bnew=zeros(npart,NS);
    part=zeros(npart,NS);
    part(:,1)=randn(npart,1);
    for j=1:npart
        xtemp=zeros(1,N*L+EX);
        xtemp(:,1:N*L)=x(i,1:N*L);
        xtemp(:,end)=1;
        xtemp(:,(N*L)+2:end-1)=b0(j,2:end);
        %draw states
        htemp=MU'+b0(j,:)*F'+part(j,:)*cQ;
        bnew(j,:)=htemp;
        xtemp(:,(N*L)+1:(N*L)+1)=htemp(:,1);
        
        %likelihood
        res1=(y(i,:)-xtemp*varcoef1);
        res2=(y(i,:)-xtemp*varcoef2);
        sigma1=iamat1*diag(Sbig.*exp(htemp(:,1)))*iamat1';
        sigma2=iamat2*diag(Sbig.*exp(htemp(:,1)))*iamat2';
        isigma1=invpd(sigma1);
        isigma2=invpd(sigma2);
        dsigma1=det(sigma1);
        dsigma2=det(sigma2);
        if e1(i)==1
            resx=res1*isigma1*res1';
            densx=(1/sqrt(dsigma1))*exp(-0.5*resx);
        else
           resx=res2*isigma2*res2';
            densx=(1/sqrt(dsigma2))*exp(-0.5*resx); 
        end
        
        
	tempdens=densx;
    if isnan(tempdens) || isinf(tempdens) || ~isreal(tempdens)
        dens(j)=exp(-100);
    else
        dens(j)=tempdens;
    end
    end
    sdens=sum(dens);
   prob=dens./sdens;
   index1=discretesample(prob,npart);
   
   b0=bnew(index1,:);
   liki=sdens/npart;
lik=lik+log(liki);
states(i,:)=sum(b0.*repmat(prob,1,NS))';
end

