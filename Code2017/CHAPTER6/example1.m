clear

addpath('gensys','functions');  %Sims (2002) model solution

%some generated data for simple 3 equation model (see pdf file)
sigma = 1;  
beta = 0.99; 
delta = 1.5;  
alpha = 3;   
omega = 1.5;  
rho1 = 0.7;     
rho2 = 0.7;     
rho3 = 0.7;  
Thetat=[sigma delta alpha omega rho1 rho2 rho3];

[  F, g, PROBLEM,GAM0,GAM1,PSI,PPI ] = model_solve( Thetat );


%%%generate Artificial Data to be used in the next example%%%%%%
t1=200;
t0=100;
t=t1+t0;
Sigma=eye(3);
ee=randn(t,3)*chol(Sigma);
ehat=zeros(t,6);
ymat=zeros(t,6);
for i=2:t
    %iid shocks
    ehat(i,:)=(g*ee(i,:)')';
    %data
    ymat(i,:)=(F*ymat(i-1,:)')'+ehat(i,:);
end

%%%%%%%%%test%%%%%%%%%%%%%%%%%%
y=ymat(t0+1:t,1:3);  %data on Y p R
Thetat=[Thetat diag(Sigma)'];
save('data','y','Thetat');





%%%Calculate IRF%%%%%%
t=20;
Sigma=eye(3);
ee=zeros(t,3);
ee(2,3)=1;
ee=ee*chol(Sigma);
ehat=zeros(t,6);
ymat=zeros(t,6);
for i=2:t
    %iid shocks
    ehat(i,:)=(g*ee(i,:)')';
    %data
    ymat(i,:)=(F*ymat(i-1,:)')'+[ehat(i,:)];
end
IRFT=ymat;
save('IRFT','IRFT')
