function out = iwpq(v,ixpx)

k=rows(ixpx);
z=zeros(v,k);
cixpx=cholx(ixpx)';
for i=1:v
    z(i,:)=(cixpx*randn(k,1))';
end
out=inv(z'*z);

