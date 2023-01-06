function [iamat1,amat1,iamat2,amat2]=getamatonlySTAR(res,Sbigin,hlast,c0,p00_a,vs0,sx0,N,e2);
amatx1=[];
amatx2=[];

for j=2:N

%v2=-a1*v1+sqrt(exp(h2))*e2
ytemp=res(:,j)+res(:,j).*e2;  %v2
Nx=cols(res(:,1:j-1));
xtemp=[res(:,1:j-1)*(-1)  (res(:,1:j-1)*(-1).*repmat(e2,1,Nx))]; %-v1
ytemp=ytemp./sqrt(hlast(1:rows(hlast),1).*Sbigin(j));
xtemp=xtemp./repmat(sqrt(hlast(1:rows(hlast),1).*Sbigin(j)),1,cols(xtemp));
%prior means and variances
alpha21_0=[c0(j,1:j-1)';c0(j,1:j-1)'];
v00_a=diag(abs(alpha21_0))*p00_a;

alpha21_2=getreg(ytemp,xtemp,alpha21_0,v00_a,1); %draw from Normal
amatx1=[amatx1 alpha21_2(1:Nx)'];
amatx2=[amatx2 alpha21_2(Nx+1:end)'];
end
amat1=repmat(amatx1,size(res,1),1);
iamat1=invpd(chofac(N,amatx1'));%chofac converts to a lower triangular 

amat2=repmat(amatx2,size(res,1),1);
iamat2=invpd(chofac(N,amatx2'));%chofac converts to a lower triangular 