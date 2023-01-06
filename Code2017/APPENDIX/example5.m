clear %clears all variables in memory
%read in data from an excel file
[data,names]=xlsread('\data\data.xls');
%dimensions of the data
T=size(data,1); %number of rows in the data
%assign data
Y=data(:,1); %Y is the first column of data
%X matrix with a constant
X=[ones(T,1) data(:,2)];
%OLS coefficients
B=inv(X'*X)*(X'*Y);
%Residuals
E=Y-X*B;
%Variance of Error term
S=(E'*E)/(T-2);
%Covariance of B
V=S.*inv(X'*X);
%standard errors
SE=diag(V).^0.5;