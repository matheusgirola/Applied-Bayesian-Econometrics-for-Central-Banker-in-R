function hnew=getsvolxx(hlast,g,mubar,sigmabar,errors,alpha,delta)

T=rows(errors);
hnew=zeros(T+1,1);

i=1;
%time period 0

hlead=hlast(i+1);
ss = sigmabar*g/(g + (delta^2)*sigmabar);   %variance
mu = ss*(mubar/sigmabar + delta*(log(hlead) - alpha)/g);







%draw from lognormal  using mu and ss
h = exp(mu + (ss^.5)*randn(1,1));
hnew(i)=h;


%time period 1 to t-1

for i=2:T
    hlead=hlast(i+1);
    hlag=hnew(i-1);
    yt=errors(i-1);  %note change
   
%mean and variance of the proposal log normal density
mu = alpha*(1-delta) + delta*(log(hlead)+log(hlag))/(1+delta^2);
ss = (g)/(1+delta^2);

%candidate draw from lognormal
htrial = exp(mu + (ss^.5)*randn(1,1));
if isnan(htrial) || htrial>1000
    accept=0;
else
%acceptance probability in logs
lp1 = -0.5*log(htrial) - (yt^2)/(2*htrial);  %numerator
lp0 = -0.5*log(hlast(i)) - (yt^2)/(2*hlast(i));   %denominator
accept = min([1;exp(lp1 - lp0)]);  %ensure accept<=1
end
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
mu = alpha + delta*log(hlag);
ss = g;
%candidate draw from lognormal
htrial = exp(mu + (ss^.5)*randn(1,1));
if isnan(htrial) || htrial>1000
    accept=0;
else
%acceptance probability
lp1 = -0.5*log(htrial) - (yt^2)/(2*htrial);
lp0 = -0.5*log(hlast(i)) - (yt^2)/(2*hlast(i));
accept = min([1;exp(lp1 - lp0)]);  %ensure accept<=1
end

u = rand(1,1);
if u <= accept;
   h = htrial;
else
   h = hlast(i);
end
hnew(i)=h;