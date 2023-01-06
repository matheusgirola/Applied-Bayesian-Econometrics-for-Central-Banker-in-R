function xx=CHECKR(x,e)

out=[];
n=cols(x);

for i=1:n
x1=x(i);
x2=e(i);
out1=0;
if x2 ~= 0
if x1<0 && x2<0
out1=1;
end;
if x1>0 && x2>0
out1=1;
end
end
out=[out out1];
end



xx=sum(out,2)==sum(abs(e),2);
