function [beta,sigma,omega,e]=getols(y,x)
N=cols(y);
K=cols(x);
T=rows(x);
beta=zeros(K,N);
e=zeros(T,N);
sigma=zeros(N,1);
omega=zeros(K,K,N);
for j=1:N
beta(:,j)=x\y(:,j);
e(:,j)=y(:,j)-x*beta(:,j);
sigma(j)=(e(:,j)'*e(:,j))/(T-K);

omega(:,:,j)=sigma(j)*invpd(x'*x);
end