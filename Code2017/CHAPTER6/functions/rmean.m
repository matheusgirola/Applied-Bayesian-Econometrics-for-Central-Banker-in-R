function out=rmean(y,horz)

NN=horz;
m1=y(1,:);
N=size(y,1);
for i=2:NN:N
%     m1=[m1 mean(y(i:i+NN-1))];
%     m2=[m2 mean(y1(i:i+NN-1))];
    m1=[m1;mean(y(1:i,:))];
    
end
out=m1