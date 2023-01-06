clear
addpath('functions');

% generate artificial data
T=200;
Y=zeros(T,1);
b0=[0.1 1.1 -0.9];
r0=0.8;
e=randn(T,1);
v=zeros(T,1);
for i=3:T
    v(i)=v(i-1)*r0+e(i);
    Y(i)=b0(1)+[Y(i-1) Y(i-2)]*b0(2:end)'+v(i);
        end
       



X=[ones(T,1) lag0(Y,1) lag0(Y,2)];

%remove missing obs
Y=Y(3:end);
X=X(3:end,:);
T=rows(X);





%step 1 set priors and starting values
%priors for B
B0=[0;0;0];
Sigma0=eye(3);
%priors for sigma2
T0=1;
D0=0.1;
%priors for rho
rho0=0;
Sigma0r=1;

%starting values
B=B0;
rho=rho0;
sigma2=1;

%parameters for the Geweke diagnostic
Pa=0.1;
Pb=0.5;

    

reps=125000;
burn=124000;
out1=[]; %will save the inflation forecast
out2=[];
out3=[];
for i=1:reps

%step 2 Sample B conditional on sigma N(M*,V*)
%remove serial correlation
ystar=Y-lag0(Y,1)*rho;
xstar=X-lag0(X,1)*rho;
ystar=ystar(2:end);
xstar=xstar(2:end,:);
M=inv(inv(Sigma0)+(1/sigma2)*(xstar'*xstar))*(inv(Sigma0)*B0+(1/sigma2)*xstar'*ystar);
V=inv(inv(Sigma0)+(1/sigma2)*(xstar'*xstar));
chck=-1;
while chck<0
B=M+(randn(1,3)*chol(V))';
b=[B(2) B(3);1   0];
ee=max(abs(eig(b)));
if ee<=1
    chck=1;
end
end

%step 3 compute rho
y=Y-X*B;
x=lag0(y,1);
y=y(2:end);
x=x(2:end);
MM=inv(inv(Sigma0r)+(1/sigma2)*(x'*x))*(inv(Sigma0r)*rho0+(1/sigma2)*x'*y);
VV=inv(inv(Sigma0r)+(1/sigma2)*(x'*x));
%draw rho but again ensure stationarity
chck=-1;
while chck<0
rho=MM+(randn(1,1)*chol(VV))';

ee=abs(rho);
if ee<=1
    chck=1;
end
end



%step 3 sample sigma2 conditional on B from IG(T1,D1);
%compute residuals
resids=ystar-xstar*B;
%compute posterior df and scale matrix
T1=T0+T;
D1=D0+resids'*resids;
%draw from IG
z0=randn(T1,1);
z0z0=z0'*z0;
sigma2=D1/z0z0;

if i>burn
    
    out1=[out1;sigma2];
    out2=[out2;B'];
    out3=[out3;rho];end
end


%compute Geweke diagnostics for convergence
Na=ceil(Pa*rows(out1));  %first sample
Nb=ceil(Pb*rows(out1));  %second sample
Nstar=rows(out1)-Nb;

%calculation for sigma
Ma=mean(out1(1:Na)); %mean of first sample
Sa=spectrum(out1(1:Na),0); %variance of first sub-sample
Mb=mean(out1(Nstar+1:end)); %mean of second sample
Sb=spectrum(out1(Nstar+1:end),0); %variance of second sub-sample
CD1=((Ma-Mb))/sqrt((Sa/Na)+(Sb/Nb));
RNE1=(std(out1)^2)/spectrum(out1,0);

%calculation for elements of coefficient matrix
CD2=[];
RNE2=[];
for j=1:cols(out2);
Ma=mean(out2(1:Na,j)); %mean of first sample
Sa=spectrum(out2(1:Na,j),0); %variance of first sub-sample
Mb=mean(out2(Nstar+1:end,j)); %mean of second sample
Sb=spectrum(out2(Nstar+1:end,j),0); %variance of second sub-sample
CD=((Ma-Mb))/sqrt((Sa/Na)+(Sb/Nb));
CD2=[CD2;CD];
RNE2=[RNE2;(std(out2(:,j))^2)/spectrum(out2(:,j),0)];
end



%calculation for rho
Ma=mean(out3(1:Na)); %mean of first sample
Sa=spectrum(out3(1:Na),0); %variance of first sub-sample
Mb=mean(out3(Nstar+1:end)); %mean of second sample
Sb=spectrum(out3(Nstar+1:end),0); %variance of second sub-sample
CD3=((Ma-Mb))/sqrt((Sa/Na)+(Sb/Nb));
RNE3=(std(out3)^2)/spectrum(out3,0);


