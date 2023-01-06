function [ Y,X,Ystar ] = prepare( data,L,tard,tarvar )
Y=data;



%take lags
X=[];
for j=1:L
X=[X lag0(data,j) ];
end
X=[X ones(rows(X),1)];

%compute threshold variable
Ystar=lag0(Y(:,tarvar),tard);


Y=Y(max([L,tard])+1:end,:);
X=X(max([L,tard])+1:end,:);
Ystar=Ystar(max([L,tard])+1:end,:);

end

