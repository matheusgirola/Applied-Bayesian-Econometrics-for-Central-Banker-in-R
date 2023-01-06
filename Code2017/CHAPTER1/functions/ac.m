function out=ac(data)
y=data;
x=lag0(data,1);
x=[x ones(rows(data),1)];
y=y(2:end);
x=x(2:end,:);
b=x\y;
out=b(1);