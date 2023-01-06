function [ h0] = getsimplevol( Y )
N=cols(Y);
T=rows(Y);
outh=zeros(T+1,N);
for j=1:N
    tmp=(diff(Y(:,j)).^2)+0.001;
    outh(:,j)=[tmp(1:2,1);tmp];
end

h01=outh;
h0=outh;
errors0=diff(Y)+0.001;
errors0=[errors0(1,:);errors0];
for i=1:50
    for j=1:N
    h01(:,j)=getsvol(h0(:,j),0.1,log(outh(1,j)),10,errors0(:,j));
    end
    h0=h01;
end

end

