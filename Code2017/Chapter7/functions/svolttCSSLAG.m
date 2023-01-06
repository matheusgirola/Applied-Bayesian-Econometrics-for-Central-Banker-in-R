function [hnew, naccept]=svolttCSSLAG(hnewlag,hleadx,hold,FF,QQ,iQQ,cQQ,varcoef,iamat,y,x,mumat,Smat)

XX=log(hnewlag);
NS=cols(hold);
N=1;
hlead=hleadx;
BINV2=iQQ+FF'*iQQ*FF;
b2=XX*FF'*iQQ+log(hlead)*iQQ*FF;
VV=invswpx(BINV2,N);
constant=((eye(cols(FF))-FF')*mumat)'*iQQ;
MM=VV*(constant+b2)';
cVV=sqrt(VV);
htrial=(exp(MM+(randn(1,NS)*cVV)'))';
htrial(N+1:NS)=exp(XX(1:NS-N)); %lagged states
if any(isnan(htrial))||any(isinf(htrial))||~isreal(htrial);
	accept=0;
else
xmat=[x(:,1:cols(x)-2) log(htrial(:,2:end)) x(:,cols(x)-1:end)];
xmatold=[x(:,1:cols(x)-2) log(hold(:,2:end)) x(:,cols(x)-1:end)];
probnew=probx5SS(y,xmat,varcoef,iamat,htrial,Smat);
probold=probx5SS(y,xmatold,varcoef,iamat,hold,Smat);
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



