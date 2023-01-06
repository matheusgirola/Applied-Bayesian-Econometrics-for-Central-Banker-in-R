function [ out ] = irfvol( beta2,f,q,N,L,EX,LH,horz )
varcoef=reshape(beta2,N*L+EX,N);
yhat=zeros(horz+L,N);
hhat=zeros(horz+L,1);
ehat=zeros(horz+L,1);
ehat(L+1)=1;

for j=L+1:horz+L
    hhat(j,:)=hhat(j-1,:)*f+ehat(j)*sqrt(q);
    xhat=[];
    for i=1:L
        xhat=[xhat yhat(j-i,:)];
    end
    for i=0:LH
        xhat=[xhat hhat(j-i,:)];
    end
    xhat=[xhat 0];
 
    yhat(j,:)=xhat*varcoef;
end
out=yhat(L+1:end,:);
