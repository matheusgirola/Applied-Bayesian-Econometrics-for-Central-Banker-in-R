function outm=ARSVOLX(data,L)
REPS=5000;
BURN=4000;
Y=data;
X=[];
for j=1:L
    X=[X lag0(data,j)];
end
X=[X ones(rows(X),1)];
Y=Y(L+1:end,:);
X=X(L+1:end,:);
T=rows(X);

N=cols(Y);
%set priors for SVOL
sigma0=std(diff(Y)).^2;

MU0=log(diag(sigma0)); %prior mean
SS0=10;
g0=0.1;  %prior scale parameter for inverse gamma
Tg0=1;     %prior degrees of freedom
g=g0;
alpha=0;
delta=1;
B0=[0.9;0];
Sigma0=eye(2).*0.01;

hlast=(diff(Y).^2)+0.0001;   %rough guess for stochastic volatility 
hlast=[hlast(1:2,:);hlast];  %rough intial guess for svol
b0=X\Y;
ss0=eye(cols(b0)).*0.1;
e0=Y-X*b0;
beta20=b0;
epsilon=e0; %starting values for VAR residuals
%50 iterations to smooth initial estimates of stochastic volatility

for m=1:50
  hnew=zeros(T+1,N);
  for i=1:N
      htemp=getsvolxx(hlast(:,i),g0,MU0(i),SS0,epsilon(:,i),0,1);
      hnew(:,i)=htemp;
      
  end
hlast=hnew;
max(hlast)
end
out=zeros(REPS-BURN,T,N);
jgibbs=1;
for igibbs=1:REPS
    
    %sample AR
    ystar=Y./sqrt(hlast(2:end));
    xstar= X./repmat(sqrt(hlast(2:end)),1,L+1);
    [beta2,problem]=getARx(ystar,xstar,b0,ss0,1);
    if problem
        beta2=beta20;
    else
        beta20=beta2;
    end
    epsilon=Y-X*beta2; %heteroscedastic
 
   
    
   
    %sample VOL
    htemp=getsvolxx(hlast,g,MU0,SS0,epsilon,alpha,delta); 
    hlast=htemp;
    
    gerrors=log(hlast)-alpha-lag0(log(hlast),1)*delta;
    g=IG(Tg0,g0,gerrors(2:end));
    yy=log(hlast);
    xx=[lag0(yy,1) ones(rows(yy),1)];
    yy=yy(2:end,:);
    xx=xx(2:end,:);
    bdraw=getARxx(yy,xx,B0,Sigma0,g);
    delta=bdraw(1);
    alpha=bdraw(2);
    
 num2str([igibbs  alpha delta] )
if igibbs>=BURN
    for t=1:T
        tmp=volatility( beta2,hlast(t+1),1,L );
        out(jgibbs,t,:)=tmp;
        
    end
    jgibbs=jgibbs+1;
end
end

outm=squeeze(median(out));