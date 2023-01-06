function [out,beta_tt,lik]=...
    getfactorwSVOL1test(data,H,F,Mux,beta0,P00,elast,hlastc, hlastw,iamatc,iamatw,Sbigc,Sbigw,...
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
R=diag(elast(i+1,:));
Q=zeros(ns,ns);
Q(1:nfact,1:nfact)=iamatw*diag(hlastw(i+1,1).*Sbigw)*iamatw';
j=nfact+1;
for jj=1:NC
    Q(j:j+nfact-1,j:j+nfact-1)=squeeze(iamatc(jj,:,:))*...
        diag(hlastc(i+1,jj).*Sbigc(jj,:))*squeeze(iamatc(jj,:,:))';
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
K=(p10*x')*ifeta;
beta11=(beta10'+K*eta')';
p11=p10-K*(x*p10);
beta_tt(i,:)=beta11;
ptt(i,:,:)=p11;
liki=-0.5*log(2*pi)-0.5*log(det(feta))+(-0.5*(eta)*inv(feta)*(eta'));
if isreal(liki) && (1-isinf(liki))
    lik=lik+liki;
else
    lik=lik-10;
end
for i=2:T
    %
  
R=diag(elast(i+1,:));
Q=zeros(ns,ns);
Q(1:nfact,1:nfact)=iamatw*diag(hlastw(i+1,1).*Sbigw)*iamatw';
j=nfact+1;
for jj=1:NC
    Q(j:j+nfact-1,j:j+nfact-1)=squeeze(iamatc(jj,:,:))*...
        diag(hlastc(i+1,jj).*Sbigc(jj,:))*squeeze(iamatc(jj,:,:))';
    j=j+nfact;
end

    %Prediction
beta10=MUx+beta11*F';
p10=F*p11*F'+Q;
yhat=(x*(beta10)')';                                               
eta=data(i,:)-yhat;
feta=(x*p10*x')+R;
ifeta=invpd(feta);
% ifeta=woodbury(iR,x,invpd(p10));
% sum(sum(ifeta-invpd(feta)))
%updating
K=(p10*x')*ifeta;
beta11=(beta10'+K*eta')';
p11=p10-K*(x*p10);
ptt(i,:,:)=p11;
beta_tt(i,:)=beta11;
liki=-0.5*log(2*pi)-0.5*logdet(feta)+(-0.5*(eta)*inv(feta)*(eta'));
if isreal(liki) && (1-isinf(liki))
    lik=lik+liki;
else
    lik=lik-10;
end
num2str(liki)
end
% Backward recursion to calculate the mean and variance of the distribution of the state
%vector
beta2 = zeros(T,ns);   %this will hold the draw of the state variable
jv=1:(NC+1)*nfact; %index of state variables to extract
wa=randn(T,ns);

i=T;  %period t
p00=squeeze(ptt(i,jv,jv)); 
beta2(i,jv)=beta_tt(i:i,jv)+(wa(i:i,jv)*cholx(p00));   %draw for beta in period t from N(beta_tt,ptt)
%periods t-1..to .1
for i=T-1:-1:1
    

Q=zeros(ns,ns);
Q(1:nfact,1:nfact)=iamatw*diag(hlastw(i+2,1).*Sbigw)*iamatw';
j=nfact+1;
for jj=1:NC
    Q(j:j+nfact-1,j:j+nfact-1)=squeeze(iamatc(jj,:,:))*...
        diag(hlastc(i+2,jj).*Sbigc(jj,:))*squeeze(iamatc(jj,:,:))';
    j=j+nfact;
end


q=Q(jv,jv);
mu=MUx(jv);
f=F(jv,:);
pt=squeeze(ptt(i,:,:));
bm=beta_tt(i:i,:)+(pt*f'*invpd(f*pt*f'+q)*(beta2(i+1:i+1,jv)-mu-beta_tt(i,:)*f')')';  
pm=pt-pt*f'*invpd(f*pt*f'+q)*f*pt;  
beta2(i:i,jv)=bm(jv)+(wa(i:i,jv)*cholx(pm(jv,jv)));  

end

out=beta2(:,jv);