clear

T=1000; %Simulate for a 1000 periods
Y=zeros(T,1);
V=randn(T,1);
RHO=0.99;
for i=2:T
temp=Y(i-1)*RHO+V(i,1);
if temp>=0
%if temp>0 && temp<1
% if temp>0 || temp>1
    Y(i)=temp;
end
end

plot(Y);
