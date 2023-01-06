function lik=loglikelihood(theta,Y,X)
%size of Time series
T=size(Y,1);
%extract parameters
beta=theta(1:2); %coefficients
sigma=theta(3)^2; %variance of the error term
E=Y-X*beta; %calculate residuals
lik=(-T/2)*log(2*pi*sigma)-(0.5*(E'*E)/sigma);
if isnan(lik) || isinf(lik) || 1-isreal(lik)
    lik=100000;
else
    lik=-lik;
end
