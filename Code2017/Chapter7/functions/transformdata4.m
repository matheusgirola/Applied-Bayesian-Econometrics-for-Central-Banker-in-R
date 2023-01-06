function outx=transformdata4(data,id,sa)
[T,N]=size(data);
out=zeros(T,N);

%SA
dataSA=getSA(data,sa);
for j=1:N
    y=dataSA(:,j);
    if id(j)==1;
    yt=pchy(y);
    elseif id(j)==4
    yt=pchy(y);
    else
    yt=y;
    end
 out(:,j)=yt;
end
%move RGDP,CPI,STI,STP to end
outx=[out(:,5:end) out(:,1:4)];