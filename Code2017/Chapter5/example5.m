%a time-varying parameter model with stochastic volatility model
clear
addpath('functions');
%Load inflation data
Y=xlsread('\data\inflation.xlsx');
Y=((log(Y)-log(lag0(Y,4))))*100;
Y=Y(5:end,:);
T=rows(Y); 
TT0=10;  %training sample

X=[lag0(Y,1) ones(T,1)];
Y=Y(2:end);
X=X(2:end,:);



%Independence metropolis hastings algorithm for svol model

%step 1 priors for g~IG(V0,T0) and initial conditions for the stochastic
%volatility
V0=0.01;  %prior scale
T0=1;       %prior degrees of freedom

Y0=Y(1:TT0);
X0=X(1:TT0,:);
B0=X0\Y0;
E0=Y0-X0*B0;
S0=(E0'*E0)/T0;
VV0=S0*inv(X0'*X0);

mubar=log(std(E0)^2);
sigmabar=10;


%step 2 set starting values for time varying coefficient beta
beta0=B0;   %state variable  b[t-1/t-1]
p00=VV0;          %variance of state variable p[t-1/t-1]

%step 3 set prior for Q
Q0=(VV0*T0)*1e-4;
Q=Q0;  %intial values



%remove training sample
Y=Y(TT0+1:end,:);
X=X(TT0+1:end,:);
T=rows(Y); 

%step 4 starting values for stochastic volatility and 
hlast=diff(Y).^2;
hlast=[hlast(1:2);hlast]+0.0001;  %small number added to ensure no value is zero
errors=diff(Y);
errors=[errors(1);errors];  %rough estimate for the errors of observation equation



g=1;
REPS=50000;
BURN=45000;
out=[];
out1=[];
out2=[];
for j=1:REPS
%step 5 data by date metropolis hastings algorithm to draw the stochastic
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
    yt=errors(i-1);  %note change
   
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
yt=errors(i-1);   
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

%step 6 draw g from the inverse Gamma distribution

gerrors=diff(log(hnew));
g=IG(T0,V0,gerrors);  %draw from the inverse Gamma distribution


%step 7 update vale of H
hlast=hnew;

%step 8 draw the time varying coefficients using CARTER and KOHN algorithm
[beta,errors]=carterkohn1(beta0',p00,hlast,Q,Y,X);

%step 9 draw Q
errorsq=diff(beta);
scaleQ=(errorsq'*errorsq)+Q0;
Q=iwpQ(T+TT0,inv(scaleQ));
%save
if j>BURN
out=[out hlast];
out1=[out1 beta(:,1)];
out2=[out2 beta(:,2)];
end
end
TT=1917.75:0.25:2011;

subplot(2,2,1);
plot(TT,[prctile(out(2:end,:)',[50 18 84])'])
legend('Estimated posterior median','lower bound','upper bound');
title('Stochastic Volatility');
axis tight

subplot(2,2,2);
plot(TT,[prctile(out1(1:end,:)',[50 18 84])'])
legend('Estimated posterior median','lower bound','upper bound');
title('Time-Varying AR(1) Coefficient');
axis tight
subplot(2,2,3);
plot(TT,[prctile(out2(1:end,:)',[50 18 84])' ])
legend('Estimated posterior median','lower bound','upper bound');
title('Time-Varying constant');
axis tight

subplot(2,2,4);
plot(TT,[prctile((out2(1:end,:)./(1-out1(1:end,:)))',[50 18 84])' Y ])
legend('Estimated posterior median','lower bound','upper bound');
title('Long Run Mean of Inflation c_{t}/(1-b{t})');
axis tight