function plot12( t,x )
plot(t,x(:,1),'*b');
hold on
plot(t,x(:,2),'--gs');
hold on
plot(t,x(:,3),'c+');
hold on
plot(t,x(:,4),'m:');
hold on
plot(t,x(:,5),'r');
hold on
plot(t,x(:,6),'k*');
hold on
plot(t,x(:,7),'y:^');
hold on
plot(t,x(:,8),'ko');
hold on
plot(t,x(:,9),'rd');
hold on
plot(t,x(:,10),'bs');
hold on
plot(t,x(:,11),'g>');
hold on
mm=nanmean(x,2);
plot(t,mm,'k','LineWidth',2);


end

