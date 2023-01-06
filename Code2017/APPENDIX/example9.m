clear

T=1000; %Simulate for a 1000 periods
Y=zeros(T,1);
V=randn(T,1);
RHO=0.99;
i=2;
while i<T
Y(i)=Y(i-1)*RHO+V(i,1);
i=i+1;
end

plot(Y);
