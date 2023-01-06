function y = getar( rho,t )
e=randn(t,1);
y=randn(t,1);
for j=2:t
    y(j)=y(j-1)*rho+e(j);
end


end

