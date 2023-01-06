function lik=...
    predictionstep(data,H,F,Mux,beta0,P00,...
    elast,hlastc, hlastw,iamatc,iamatw,Sbigc,Sbigw,...
    NC,nfact)
T=rows(data);
ns=cols(P00);

%Carter and Kohn algorithm to draw the factor
beta_tt=zeros(T,ns);          %will hold the filtered state variable
ptt=zeros(T,ns,ns);    % will hold its variance
% %%%%%%%%%%%Step 6a run Kalman Filter
lik=0;
i=1;
x=H;
MUx=Mux';
R=diag(exp(elast));
Q=zeros(ns,ns);
Q(1:nfact,1:nfact)=iamatw*diag(exp(hlastw).*Sbigw)*iamatw';
j=nfact+1;
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
ifeta=invpd(feta);
% ifeta=woodbury(iR,x,invpd(p10));
% sum(sum(ifeta-invpd(feta)))
%updating
% K=(p10*x')*ifeta;
% beta11=(beta10'+K*eta')';
% p11=p10-K*(x*p10);
% beta_tt(i,:)=beta11;
% ptt(i,:,:)=p11;
liki=-0.5*log(2*pi)-0.5*log(det(feta))+(-0.5*(eta)*ifeta*(eta'));
if isreal(liki) && (1-isinf(liki))
    lik=liki;
else
    lik=-10;
end
