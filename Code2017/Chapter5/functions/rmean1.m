function out=rmean(y,horz)

NN=horz;
m1=[];
N=size(y,1);
for i=1:NN:N-NN
temp=y(i:i+NN-1,:,:);
temp1=mean(temp,1);
temp1=squeeze(temp1);

    m1=[m1;vec(temp1)'];
    
end
out=m1;