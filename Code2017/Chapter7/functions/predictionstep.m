function lik=...
    predictionstep(data,fload,F,Mux,beta0,P00,...
    elast,hlastc, hlastw,hlaste,iamatc,iamatw,iamate,Sbigc,Sbigw,Sbige,...
    NC,nfact,rho,index,id,L)
T=rows(data);
ns=cols(P00);
NN=cols(data);

% %%%%%%%%%%%Step 6a run Kalman Filter
lik=0;
i=1;
x=getH( fload,nfact,NN,NC,index,id,rho,L );

MUx=Mux';
R=diag(exp(elast));
Q=zeros(ns,ns);
Q(1:nfact,1:nfact)=iamatw*diag(exp(hlastw).*Sbigw)*iamatw';
Q(nfact+1:nfact*2,nfact+1:nfact*2)=iamate*diag(exp(hlaste).*Sbige)*iamate';

j=(nfact*2)+1;
for jj=1:NC
    Q(j:j+nfact-1,j:j+nfact-1)=squeeze(iamatc(jj,:,:))*...
        diag(exp(hlastc(jj)).*Sbigc(jj,:))*squeeze(iamatc(jj,:,:))';
    j=j+nfact;
end

%Prediction
beta10=MUx+beta0*F';
p10=F*P00*F'+Q;
yhat=(x*(beta10)')';                                                
eta=data(i,:)-yhat;

feta=(x*p10*x')+R;
ifeta=pinv(feta);%invpd(feta);

liki=-0.5*log(2*pi)-0.5*log(det(feta))+(-0.5*(eta)*ifeta*(eta'));
if isreal(liki) && (1-isinf(liki))
    lik=liki;
else
    lik=-10;
end
