function [ Yw,Xw ] = preparex( dataw,L,CC )
Yw=packr(dataw);
Xw=[];
for j=1:L
Xw=[Xw lag0(dataw,j) ];
end
if CC==1
Xw=[Xw ones(rows(Xw),1) ];
end

Yw=Yw(L+1:end,:);
Xw=Xw(L+1:end,:);

end

