function outm=AR0SVOL(data)
REPS=5000;
BURN=4000;
Y=data;

T=rows(Y);

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

epsilon=Y; %starting values for VAR residuals
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
    
    
    
    %sample VOL
    htemp=getsvol(hlast,g,MU0,SS0,Y); 
    hlast=htemp;
    
    gerrors=diff(log(hlast));
    g=IG(Tg0,g0,gerrors);
 igibbs   
if igibbs>=BURN
    out(jgibbs,:,:)=hlast(2:end);
    jgibbs=jgibbs+1;
end
end

outm=squeeze(median(out));