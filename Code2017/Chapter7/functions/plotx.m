function h=plotx(x)
h=plot(x(:,1),'r','LineWidth',2)
hold on
plot(x(:,2:3),'k.-');
grid on
axis tight