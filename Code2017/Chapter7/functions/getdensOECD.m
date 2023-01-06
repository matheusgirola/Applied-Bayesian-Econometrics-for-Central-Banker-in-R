
function[dens,hneww,hnewc,rnew,pmatnew]= getdensOECD(b0w,parthw,...
    b0c,parthc,b0r,partr,...
    Sbigw,Sbigc,b0p,partp,y,MUw,Fw,MUc,Fc,H,MUr,Fr,Qr,...
    nfact,Fx,Mux,iAAw,iAAc,NS,NC)
%draw states

        rtemp=MUr'+b0r.*Fr'+partr.*sqrt(Qr');
        
        rnew=rtemp;
        
        htempw=MUw'+b0w*Fw'+parthw;
        hneww=htempw;

        htempc=MUc'+b0c*Fc'+parthc;

        hnewc=htempc;
        %volatilities
        sigmaw=iAAw*diag(Sbigw.*exp(htempw(1)))*iAAw';
        csigmaw=cholx(sigmaw);
        Q=zeros(NS,NS);
Q(1:nfact,1:nfact)=csigmaw;
% jj=1;
% squeeze(iAAc(jj,:,:))
j=nfact+1;
for jj=1:NC
    Q(j:j+nfact-1,j:j+nfact-1)=cholx(squeeze(iAAc(jj,:,:))*...
        diag(exp(htempc(jj)).*Sbigc(jj,:))*squeeze(iAAc(jj,:,:))');
    j=j+nfact;
end
        
     
    
%         size(partp)
%         size(Q)

        pmattemp=Mux'+b0p*Fx'+partp*Q;
        
        pmatnew=pmattemp;
        
        %likelihood
       
        yhat=(pmattemp*H');
       
        res=y-yhat;
        RR=diag(exp(rtemp));
        dR=logdet(RR);
        iRR=diag(1./exp(rtemp));
        
            resx=res*iRR*res';
%             densx=(1/sqrt(dR))*exp(-0.5*resx);
denslog=-0.5*dR+(-0.5*resx)    ;

tempdens=(denslog);
    if isnan(tempdens) || isinf(tempdens) || ~isreal(tempdens)
        dens=exp(-100);
    else
        dens=tempdens;
    end