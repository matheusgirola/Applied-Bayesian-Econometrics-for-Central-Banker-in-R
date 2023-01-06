function plotx8(t,y)
set(gcf,'DefaultAxesColorOrder',[0 0 1; 1 0 0;1 0 0;0 0 1]);
cu=y(:,2);
cl=y(:,1);

h=t;
h=h';
hh=fill([h(1); h(1:end); flipud([h(1:end); h(end)])],[cu(1); cl(1:end); flipud([cu(1:end); cl(size(cl,1))])],'b');
% set(hh,'edgecolor',[0.84 0.93 1]);
% set(hh,'facecolor',[0.84 0.93 1]);
set(hh,'edgecolor',[0 0.85 0.85]);
set(hh,'facecolor',[0 0.85 0.85]);
% set(hh,'edgecolor',[0.5 0.85 0.85]);
% set(hh,'facecolor',[0.5 0.85 0.85]);

% hold on;
% 
% plot(h,y(:,1),'LineWidth',2);