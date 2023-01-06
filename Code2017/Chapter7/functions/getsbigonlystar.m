function [Sbig]=getsbigonlystar(res,Sbigin,hlast,vs0,sx0,N,A1,A2,e2);
Sbig=Sbigin;
AA1=chofac(N,A1(1,:)');
AA2=chofac(N,A2(1,:)');
e1=1-e2;
for j=2:N
%sample Sbig

ytemp=res(:,j)+res(:,j).*e2;  %v2
xtemp=res(:,1:j-1)*(-1); %-v1
ytemp=ytemp./sqrt(hlast(1:rows(hlast),1));
xtemp=xtemp./repmat(sqrt(hlast(1:rows(hlast),1)),1,cols(xtemp));
alpha21_1=AA1(j,1:j-1)';
alpha21_2=AA2(j,1:j-1)';
resx=((((ytemp-xtemp*alpha21_1).*e1)-((xtemp*alpha21_2).*e2)));
Sbig(j)=IG(vs0,sx0(j),resx);
end
