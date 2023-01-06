
function[dens,hnew,rnew,pmatnew]= getdensnorho(b0,parth,...
    b0r,partr,...
    Sbig,b0p,partp,N,L,EX,ns,y,MU,F,H,mumatr,fmatr,Qr,...
    nfact,varcoef,iAA,bload,nsr,NN)
%draw states
rnew=zeros(NN,nsr);
exor=zeros(1,NN);
RR=zeros(NN,NN);
iRR=RR;
for ii=1:NN
Fr=zeros(nsr,nsr);
Fr(1:1,1:1)=fmatr(ii);
Fr(2:nsr,1:nsr-1)=eye(nsr-1);

cQr=zeros(nsr,nsr);
cQr(1,1)=sqrt(Qr(ii));
MUr=zeros(nsr,1);
MUr(1)=mumatr(ii);
rlag=squeeze(b0r(1,ii,:))';
rtemp=MUr'+rlag*Fr'+[partr(:,ii) zeros(1,nsr-1)]*cQr;
rnew(ii,:)=rtemp;

exor(:,ii)=rtemp*bload(ii,:)'; 
RR(ii,ii)=exp(rtemp(1));
iRR(ii,ii)=1./exp(rtemp(1));
end
        
        
        htemp=MU'+b0*F'+parth;
        hnew=htemp;
        
        %volatilities
        
        sigma=iAA*diag(Sbig.*exp(htemp(:,1)))*iAA';
        csigma=zeros(ns,ns);
        csigma(1:nfact,1:nfact)=cholx(sigma);
        xtemp=zeros(1,N*L+EX);
        xtemp(:,(N*L)+2:end-1)=b0(:,2:end);
        xtemp(:,(N*L)+1:(N*L)+1)=htemp(:,1);
        xtemp(:,end)=1;
        exo=xtemp(:,N*L+1:end)*varcoef(N*L+1:end,:);
        mup=zeros(1,ns);
        mup(:,1:nfact)=exo;
        fp=comp(vec(varcoef),N,L,EX);
        
      
        pmattemp=mup+b0p*fp'+partp*csigma;
        
        pmatnew=pmattemp;
        
        %likelihood
        
        yhat=(pmattemp*H');
       
        res=y-yhat-exor;
        RR=diag(exp(rtemp));
        dR=logdet(RR);
       
        
            resx=res*iRR*res';
%             densx=(1/sqrt(dR))*exp(-0.5*resx);
denslog=-0.5*dR+(-0.5*resx)    ;

tempdens=(denslog);
    if isnan(tempdens) || isinf(tempdens) || ~isreal(tempdens)
        dens=exp(-100);
    else
        dens=tempdens;
    end