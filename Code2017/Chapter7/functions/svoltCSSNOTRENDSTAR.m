function [hnew, naccept]=svoltCSSNOTRENDSTAR(hnewlag,hold,FF,QQ,iQQ,cQQ,varcoef1,varcoef2,iamat1,iamat2,y,x,mumat,Sbig,e2)

XX=log(hnewlag);
NS=cols(hold);
N=1;
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
xmat=[x(:,1:cols(x)-1) log(htrial) x(:,end:end)];
xmatold=[x(:,1:cols(x)-1) log(hold) x(:,end:end)];
probnew=probx5SSSTAR(y,xmat,varcoef1,varcoef2,iamat1,iamat2,htrial,Sbig,e2);
probold=probx5SSSTAR(y,xmatold,varcoef1,varcoef2,iamat1,iamat2,hold,Sbig,e2);
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



