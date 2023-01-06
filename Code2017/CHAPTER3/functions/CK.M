function[beta2,lik]=ck(T,ns,MU,F,Q,R,H,datastar,beta0,P00)

%step 7 Carter and Kohn algorithm to draw the factor
%Carter and Kohn algorithm to draw the factor
beta_tt=[];          %will hold the filtered state variable
ptt=zeros(T,ns,ns);    % will hold its variance
% %%%%%%%%%%%Step 6a run Kalman Filter
lik=0;
i=1;
x=H;
%Prediction
beta10=MU+beta0*F';
p10=F*P00*F'+Q;
yhat=(x*(beta10)')';                                                
eta=datastar(i,:)-yhat;
feta=(x*p10*x')+R;
%updating
K=(p10*x')*inv(feta);
beta11=(beta10'+K*eta')';
p11=p10-K*(x*p10);
beta_tt=[beta_tt;beta11];
ptt(i,:,:)=p11;
liki=-0.5*log(2*pi)-0.5*log(det(feta))+(-0.5*(eta)*inv(feta)*(eta'));
lik=lik+liki;
for i=2:T
    %Prediction
beta10=MU+beta11*F';
p10=F*p11*F'+Q;
yhat=(x*(beta10)')';                                               
eta=datastar(i,:)-yhat;
feta=(x*p10*x')+R;
%updating
K=(p10*x')*inv(feta);
beta11=(beta10'+K*eta')';
p11=p10-K*(x*p10);
ptt(i,:,:)=p11;
beta_tt=[beta_tt;beta11];
liki=-0.5*log(2*pi)-0.5*log(det(feta))+(-0.5*(eta)*inv(feta)*(eta'));
lik=lik+liki;
end
% Backward recursion to calculate the mean and variance of the distribution of the state
%vector
beta2 = zeros(T,ns);   %this will hold the draw of the state variable
jv=1; %index of state variables to extract
wa=randn(T,ns);
f=F(jv,:);
q=Q(jv,jv);
mu=MU(jv);
i=T;  %period t
p00=squeeze(ptt(i,jv,jv)); 
beta2(i,jv)=beta_tt(i:i,jv)+(wa(i:i,jv)*cholx(p00));   %draw for beta in period t from N(beta_tt,ptt)
%periods t-1..to .1
for i=T-1:-1:1
pt=squeeze(ptt(i,:,:));
bm=beta_tt(i:i,:)+(pt*f'*inv(f*pt*f'+q)*(beta2(i+1:i+1,jv)-mu-beta_tt(i,:)*f')')';  %update the filtered beta for information contained in beta[t+1] 
                                                                               %i.e. beta2(i+1:i+1,:) eq 8.16 pp193 in Kim Nelson
pm=pt-pt*f'*inv(f*pt*f'+q)*f*pt;  %update covariance of beta
beta2(i:i,jv)=bm(jv)+(wa(i:i,jv)*cholx(pm(jv,jv)));  %draw for beta in period t from N(bm,pm)eq 8.17 pp193 in Kim Nelson

end