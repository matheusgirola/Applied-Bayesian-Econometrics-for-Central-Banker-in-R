%a stochastic volatility model for UK inflation
clear

addpath('functions');
%oad inflation data
Y=xlsread('\data\inflation.xlsx');
Y=((log(Y)-log(lag0(Y,4))))*100;
Y=Y(5:end,:);
T=rows(Y); 
TT0=10;  %training sample
    
%Independence metropolis hastings algorithm for svol model

%step 1 priors for g~IG(V0,T0) and initial conditions for the stochastic
%volatility
V0=0.01;  %prior scale
T0=1;       %prior degrees of freedom
mubar=log(std(Y(1:TT0))^2);
sigmabar=10;


%remove training sample
Y=Y(TT0+1:end,:);
T=rows(Y); 

%step 2 starting values for stochastic volatility and 
hlast=diff(Y).^2;
hlast=[hlast(1:2);hlast]+0.0001;  %small number added to ensure no value is zero

g=1;
REPS=30000;
BURN=25000;
out=[];
for j=1:REPS
%step 3 data by date metropolis hastings algorithm to draw the stochastic
%volatility
hnew=zeros(T+1,1);


i=1;
%time period 0

hlead=hlast(i+1);
ss = sigmabar*g/(g + sigmabar);   %variance
mu = ss*(mubar/sigmabar + log(hlead)/g);  %mean
%draw from lognormal  using mu and ss
h = exp(mu + (ss^.5)*randn(1,1));
hnew(i)=h;


%time period 1 to t-1

for i=2:T
    hlead=hlast(i+1);
    hlag=hnew(i-1);
    yt=Y(i-1);
   
%mean and variance of the proposal log normal density
mu = (log(hlead)+log(hlag))/2;  
ss = g/2;

%candidate draw from lognormal
htrial = exp(mu + (ss^.5)*randn(1,1));

%acceptance probability in logs
lp1 = -0.5*log(htrial) - (yt^2)/(2*htrial);  %numerator
lp0 = -0.5*log(hlast(i)) - (yt^2)/(2*hlast(i));   %denominator
accept = min([1;exp(lp1 - lp0)]);  %ensure accept<=1

u = rand(1,1);
if u <= accept;
   h = htrial;
else
   h = hlast(i);
end
hnew(i)=h;
end


%time period T
i=T+1;
yt=Y(i-1);
hlag=hnew(i-1);


%mean and variance of the proposal density
mu = log(hlag);   % only have ht-1
ss = g;
%candidate draw from lognormal
htrial = exp(mu + (ss^.5)*randn(1,1));

%acceptance probability
lp1 = -0.5*log(htrial) - (yt^2)/(2*htrial);
lp0 = -0.5*log(hlast(i)) - (yt^2)/(2*hlast(i));
accept = min([1;exp(lp1 - lp0)]);  %ensure accept<=1


u = rand(1,1);
if u <= accept;
   h = htrial;
else
   h = hlast(i);
end
hnew(i)=h;

%step 4 draw g from the inverse Gamma distribution

errors=diff(log(hnew));
g=IG(T0,V0,errors);  %draw from the inverse Gamma distribution


%step 5 update vale of H
hlast=hnew;

%save
if j>BURN
out=[out hlast];
end
end

TT=1917.5:0.25:2011;
subplot(1,2,1);
plot(TT(1:end),Y);
title('Annual CPI inflation for the UK');
axis tight
subplot(1,2,2);
plot(TT,[prctile(out(2:end,:)',[50 18 84])']);
title('Estimated stochastic volatility');
axis tight


legend('Estimated posterior median','lower bound','upper bound','true');