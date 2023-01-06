function [AX,H1,H2] = plotvolx( t,z1,t1,data,names )
[AX,H1,H2]=plotyy(t,z1,t1,data,'plotx2','plotkk');

% title(strcat('Probability of high volatility regime',':',names));
legend(names)
set(AX(1),'xlim',[min(t) max(t)]);
set(AX(2),'xlim',[min(t) max(t)]);
set(AX(1),'ylim',[min(min(z1)) max(max(z1))]);
set(AX(2),'ylim',[min(min(data)) max(max(data))]);

