function [hnew, naccept]=svoltC(hnewlag,hold,FF,QQ,iQQ,cQQ,varcoef,iamat,y,x,mumat)

XX=log(hnewlag);
NS=cols(hold);
N=cols(y);
BINV2=iQQ+FF'*iQQ*FF;
b2=XX*FF'*iQQ;
VV=invswpx(BINV2,N);
constant=((eye(cols(FF))-FF')*mumat)'*iQQ;
MM=VV*(constant+b2)';
cVV=sqrt(VV);
htrial=(exp(MM+(randn(1,NS)*cVV)'))';
htrial(N+1:NS)=exp(XX(1:NS-N)); %lagged states
if any(isnan(htrial))||any(isinf(htrial))||~isreal(htrial);
	accept=0;
else
xmat=[x(:,1:cols(x)-1) log(htrial) x(:,cols(x))];
xmatold=[x(:,1:cols(x)-1) log(hold) x(:,cols(x))];
probnew=probx5(y,xmat,varcoef,iamat,htrial);
probold=probx5(y,xmatold,varcoef,iamat,hold);
accept=exp(probnew-probold);

end
hnew=zeros(rows(hold),cols(hold));
u=rand(1,1);
if u<accept;
    hnew(:,1:N)=htrial(:,1:N);
    hnew(:,N+1:NS)=hnewlag(:,1:NS-N);
    naccept=1;
else
	  hnew(:,1:N)=hold(:,1:N);
    hnew(:,N+1:NS)=hnewlag(:,1:NS-N);
	naccept=0;
end



