clear
load output

figure(1) %plot raw Gibbs sample output

subplot(3,2,1)
plot(out2(:,1))
title('\alpha');
xlabel('Gibbs iterations')
axis tight
subplot(3,2,2)
plot(out2(:,2))
title('B_{1}');
xlabel('Gibbs iterations')
axis tight

subplot(3,2,3)
plot(out2(:,3))
title('B_{2}');
xlabel('Gibbs iterations')
axis tight
subplot(3,2,4)
plot(out3(:,1))
title('\rho');
xlabel('Gibbs iterations')
axis tight


subplot(3,2,5)
plot(out4(:,1))
title('\sigma');
xlabel('Gibbs iterations')
axis tight

figure(2) %plot recursive means of Gibbs sample output
H=20;
subplot(3,2,1)
plot(rmean1(out2(:,1),H))
title('\alpha');
xlabel('Gibbs iterations')
axis tight

subplot(3,2,2)
plot(rmean1(out2(:,2),H))
title('B_{1}');
xlabel('Gibbs iterations')
axis tight


subplot(3,2,3)
plot(rmean1(out2(:,3),H))
title('B_{2}');
xlabel('Gibbs iterations')
axis tight

subplot(3,2,4)
plot(rmean1(out3(:,1),H))
title('\rho');
xlabel('Gibbs iterations')
axis tight

subplot(3,2,5)
plot(rmean1(out4(:,1),H))
title('\sigma');
xlabel('Gibbs iterations')
axis tight

figure(3) %plot raw Gibbs sample output

subplot(3,2,1)
autocorr(out2(:,1));
title('\alpha');

axis tight
subplot(3,2,2)
autocorr(out2(:,2));
title('B_{1}');

axis tight

subplot(3,2,3)
autocorr(out2(:,3));
title('B_{2}');


axis tight
subplot(3,2,4)
autocorr(out3(:,1));
title('\rho');

axis tight


subplot(3,2,5)
autocorr(out4(:,1));
title('\sigma');
xlabel('Gibbs iterations')
axis tight


