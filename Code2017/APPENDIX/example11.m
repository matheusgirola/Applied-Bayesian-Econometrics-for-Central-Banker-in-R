
clear %clears all variables in memory

addpath('func') %adds the functions in folder func
global Y X;
%read in data from an excel file
[data,names]=xlsread('\data\data.xls');
%dimensions of the data
T=size(data,1); %number of rows in the data
%assign data
Y=data(:,1); %Y is the first column of data
%X matrix with a constant
X=[ones(T,1) data(:,2)];
%specify starting values
theta=[0;0;0.1];
%OLS for comparison
bols=X\Y;
eols=Y-X*bols;
sols=(eols'*eols)/T;
% call csminwel
[fhat,xhat,ghat,Hhat,itct,fcount,retcodehat] = csminwel('loglikelihood',theta,eye(3)*.5,[] ,1e-14,100,Y,X);

%try simulated annealing

Bounds=[ones(length(theta),1).*-10 ones(length(theta),1).*10];
soptions = saoptimset('Display','iter','Tolfun',1.0e-5,'InitialTemperature',50,'MaxIter',55000,'StallIterLimit',1000);
soptions=saoptimset(soptions,'PlotFcns',{@saplotbestf,@saplottemperature,@saplotf,@saplotstopping},'ReannealInterval', 50);

 [xopt,fval,exitflag,output] = simulannealbnd(@loglikelihoodSA,theta,Bounds(:,1),Bounds(:,2),soptions);

disp('---------------------')
disp('OLS          CSMINWEL      SA')
disp([[bols;sols] [xhat(1:2);xhat(3)^2] [xopt(1:2);xopt(3)^2]]);
disp('----------------------');