function x=vecr(y)
x=[];
for i=1:rows(y)
    x=[x;y(i,:)'];
end