function outm=ARSVOL(data,L)
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
g0=0.01^2;  %prior scale parameter for inverse gamma
Tg0=1;     %prior degrees of freedom
g=g0;
hlast=(diff(Y).^2)+0.0001;   %rough guess for stochastic volatility 
hlast=[hlast(1:2,:);hlast];  %rough intial guess for svol
b0=X\Y;

e0=Y-X*b0;
epsilon=e0; %starting values for VAR residuals
%50 iterations to smooth initial estimates of stochastic volatility

for m=1:50
  hnew=zeros(T+1,N);
  for i=1:N
      htemp=getsvol(hlast(:,i),g0,MU0(i),SS0,epsilon(:,i));
      hnew(:,i)=htemp;
      
  end
hlast=hnew;
end
out=zeros(REPS-BURN,T,N);
jgibbs=1;
for igibbs=1:REPS
    
    %sample AR
    ystar=Y./sqrt(hlast(2:end));
    xstar= X./repmat(sqrt(hlast(2:end)),1,L+1);
    M=xstar\ystar;
    V=invpd(xstar'*xstar);
    chck=-1;
    while chck<0
    beta2=M+(randn(1,L+1)*cholx(V))';
    S=stability(beta2,N,L,1);
    if S==0
        chck=1;
    end
    end
    epsilon=Y-X*beta2; %heteroscedastic
 
    
    %sample VOL
    htemp=getsvol(hlast,g,MU0,SS0,epsilon); 
    hlast=htemp;
    
    gerrors=diff(log(hlast));
    g=IG(Tg0,g0,gerrors);
 igibbs   
if igibbs>=BURN
    for t=1:T
        tmp=volatility( beta2,hlast(t+1),1,L );
        out(jgibbs,t,:)=tmp;
        
    end
    jgibbs=jgibbs+1;
end
end

outm=squeeze(median(out));