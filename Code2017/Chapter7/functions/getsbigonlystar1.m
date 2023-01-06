function [Sbig]=getsbigonlystar1(res,Sbigin,hlast,vs0,sx0,N,A);
Sbig=Sbigin;
AA1=chofac(N,A(1,:)');

for j=2:N
%sample Sbig

ytemp=res(:,j);  %v2
xtemp=res(:,1:j-1)*(-1); %-v1
ytemp=ytemp./sqrt(hlast(1:rows(hlast),1));
xtemp=xtemp./repmat(sqrt(hlast(1:rows(hlast),1)),1,cols(xtemp));
alpha21_1=AA1(j,1:j-1)';

resx=((ytemp-xtemp*alpha21_1));
Sbig(j)=IG(vs0,sx0(j),resx);
end
