function Y=AR(RHO,T)
Y=zeros(T,1);
V=randn(T,1);
for i=2:T
Y(i)=Y(i-1)*RHO+V(i,1);
end
