function[beta2,error,roots,epsilon]=cohnvar(Y,X,Q,amat,beta_tt,ptt,L,randmat)
T=rows(Y);
N=cols(Y);
%%Step 2a Set up matrices for the Kalman Filter

ns=cols(beta_tt);
F=eye(ns);
mu=0;


%%%%%%%%%%%end of Kalman Filter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%step 2c Backward recursion to calculate the mean and variance of the distribution of the state
%vector

    
beta2 = zeros(T,ns);   %this will hold the draw of the state variable
wa=randmat;
error=zeros(T,N);
epsilon=zeros(T,N);
roots=zeros(T,1);

i=T;  %period t
p00=squeeze(ptt(i,:,:)); 
beta2(i,:)=beta_tt(i:i,:)+(wa(i:i,:)*cholx(p00));   %draw for beta in period t from N(beta_tt,ptt)
error(i,:)=Y(i,:)-X(i,:)*reshape(beta2(i:i,:),N*L+1,N);  %var residuals
a=amat(i,:);
A=chofac(N,a');
epsilon(i,:)=error(i,:)*A';
roots(i)=stability(beta2(i,:)',N,L);

%periods t-1..to .1

for i=T-1:-1:1
pt=squeeze(ptt(i,:,:));
iFptF=invpd(F*pt*F'+Q);
% size(iFptF)
bm=beta_tt(i:i,:)+(pt*F'*iFptF*(beta2(i+1:i+1,:)-beta_tt(i,:)*F')')';  %update the filtered beta for information contained in beta[t+1]                                                                                 %i.e. beta2(i+1:i+1,:) eq 8.16 pp193 in Kim Nelson
pm=pt-pt*F'*iFptF*F*pt;  %update covariance of beta
beta2(i:i,:)=bm+(wa(i:i,:)*cholx(pm));  %draw for beta in period t from N(bm,pm)eq 8.17 pp193 in Kim Nelson
error(i,:)=Y(i,:)-X(i,:)*reshape(beta2(i:i,:),N*L+1,N);  %var residuals
a=amat(i,:);
A=chofac(N,a');
epsilon(i,:)=error(i,:)*A';
roots(i)=stability(beta2(i,:)',N,L);

end