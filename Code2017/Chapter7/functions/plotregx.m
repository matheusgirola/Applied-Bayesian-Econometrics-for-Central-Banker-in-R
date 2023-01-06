function [AX,H1,H2] = plotregx( t,z1,data,names )
[AX,H1,H2]=plotyy(t,z1,t,data,'bar','plot');

% title(strcat('Probability of high volatility regime',':',names));
legend(names)
set(AX(1),'xlim',[min(t) max(t)]);
set(AX(2),'xlim',[min(t) max(t)]);
set(AX(1),'ylim',[0 1]);
set(AX(2),'ylim',[min(min(data)) max(max(data))]);
set(AX(1),'Ycolor',[0.95 0.95 0.95]);
set(AX(2),'Ycolor','k');
set(H1,'FaceColor',[0.85 0.85 0.85]);
set(H1,'EdgeColor',[0.85 0.85 0.85]);
