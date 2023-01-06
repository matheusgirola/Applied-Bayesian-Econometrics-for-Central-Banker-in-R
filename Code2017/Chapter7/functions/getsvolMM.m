function[hnew,naccept]= ...
    getsvolMM(fbig,Qbig,mubig,betaf,hlast,mu0,is0,y,x)
NS=cols(is0);
T=rows(y);
naccept=zeros(T+1,1);
hnew=zeros(T+1,NS);
%%%%%%%%%%%%Step 5 Stochastic vol%%%%%%%%
Fmat=zeros(NS,NS);
Fmat(1:1,1:1)=fbig(1);
Qmat=zeros(NS,NS);
Qmat(1:1,1:1)=Qbig;
iQmat=zeros(NS,NS);
iQmat(1:1,1:1)=1/(Qbig);
mumat=zeros(NS,1);
mumat(1,1)=mubig(1,1);
Fmat(2:NS,1:NS-1)=eye(NS-1);
cQQ=sqrt(Qmat);
varcoef=betaf'; 
tt=1;                     %initial period
hleadx=hlast(tt+1,:);
[hh,nc]=svol1C(mu0,is0,hleadx,Fmat,Qmat,iQmat,mumat); %see Carlin Polson Stoffer example 1.1 expression for t=0


hnew(tt,:)=hh;
naccept(tt)=naccept(tt)+nc;

%time period 1 to T
for tt=2:T
hleadx=hlast(tt+1,:);
holdx=hlast(tt,:);

%////independence metropolis using f(h[t]\h[t+1],h[t-1]) as the candidate density
[hh,nc]=svolttCSSNOTRENDM(hnew(tt-1,:),hleadx,holdx,Fmat,Qmat,iQmat,cQQ,varcoef,y(tt-1,:),x(tt-1,:),mumat);

hnew(tt,:)=hh;
naccept(tt)=naccept(tt)+nc;
end

%time period T+1
tt=T+1;
holdx=hlast(tt,:);

[hh,nc]=svoltCSSNOTRENDM(hnew(tt-1,:),holdx,Fmat,Qmat,iQmat,cQQ,varcoef,y(tt-1,:),x(tt-1,:),mumat);

hnew(tt,:)=hh;
naccept(tt)=naccept(tt)+nc;
