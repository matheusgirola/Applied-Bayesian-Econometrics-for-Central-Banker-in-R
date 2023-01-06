function [ h0] = getsimplevolx( Y )
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
    parfor j=1:N
    h01(:,j)=getsvolx(h0(:,j),0.1,log(0.1),10,errors0(:,j),0,1);
    end
    h0=h01;
end

end

