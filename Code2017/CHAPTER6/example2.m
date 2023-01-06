clear

addpath('gensys','functions','distributions');  

load data
%append Thetat with the variances of the three shocks
theta=[Thetat';1;1;1];
loglik=likelihood(theta,y);  %calculate log likelihood
logpr=logprior(theta); %evaluate prior
logpost=loglik+logpr; %log posterior