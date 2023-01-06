function out=posterior(theta,y,bounds,neg)

%check if parameters withing bounds
lb=bounds(:,1);
ub=bounds(:,2);
check=sum(theta'<=lb| theta'>ub)>0;
%calculate log likelihood
lik=likelihood(theta,y);
%calculate log prior
prior=logprior(theta);
%calculate posterior
posterior=lik+prior;

if neg==1
%return negative of posterior (useful when using a minimiser in Matlab)
if isreal(posterior) && (~isinf(posterior))  && (~check )
    out=-posterior;
else
    out=inf;
end
else
%return posterior
if isreal(posterior) && (~isinf(posterior))  && (~check )
    out=posterior;
else
    out=-inf;
end
end

