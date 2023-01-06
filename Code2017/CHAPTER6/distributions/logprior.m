function out=logprior(theta)
% sigma = 1;  %>0  gamma
% beta = 0.99; %fixed point
% delta = 1.5;  % bigger than 1.01  % gamma
% alpha = 3;    %bigger than 0  gamma
% omega = 1.5;  %bigger than 1  gamma
% rho1 = 0;     % between 0 and 1  beta
% rho2 = 0;     %between 0 and 1   beta
% rho3 = 0.5;   %between 0 an 1    beta
% Sigma1       %greater than 0 inverse gamma
% Sigma2        %inverse gamma
% Sigma3        %inverse gamma





out=0;
%sigma has a gamma prior with amean of 1 and variance 1
%compute shape and scale parameters
        mu=1;
        v=1;
        b = v/mu;
        a = mu/b;
        outi=lpdfgam(theta(1),a,b);
        out=out+outi;
        
 %delta hashas a gamma prior with amean of 1.5 and variance 1
%compute shape and scale parameters
        mu=1.5;
        v=1;
        b = v/mu;
        a = mu/b;
        
        outi=lpdfgam(theta(2),a,b);
        out=out+outi;
         
  %alpha has gamma prior with mean 3 and variance 1
  
       mu=3;
        v=1;
        b = v/mu;
        a = mu/b;
        outi=lpdfgam(theta(3),a,b);
        out=out+outi;
  
  
  %omega is Gamma with mean 1.5 variance 1
   mu=1.5;
        v=1;
        b = v/mu;
        a = mu/b;
        
        outi=lpdfgam(theta(4),a,b);
        out=out+outi;
   
   
   %rho is beta with mean 0.5 variance 0.2
     mu=0.5;
     v=0.2;
     a = (1-mu)*mu^2/v^2 - mu;
        b = a*(1/mu - 1);
        out =out + lpdfbeta(theta(5),a,b);
        
        
        %rho is beta with mean 0.5 variance 0.2
     mu=0.5;
     v=0.2;
     a = (1-mu)*mu^2/v - mu;
        b = a*(1/mu - 1);
        out =out + lpdfbeta(theta(6),a,b);
        
        
        %rho is beta with mean 0.5 variance 0.2
     mu=0.5;
     v=0.2;
     a = (1-mu)*mu^2/v - mu;
        b = a*(1/mu - 1);
        out =out + lpdfbeta(theta(7),a,b);
        
        
        %sigma is inverse Gamma(1,0.5)
        mu=1;
        v=0.5;
        mu2=mu^2;
          nu   = 2*(2+mu2/v);
    s    = 2*mu*(1+mu2/v);
        out=out+lpdfig2(theta(8),s,nu);
  
   %sigma is inverse Gamma(1,0.5)
       mu=1;
        v=0.5;
        mu2=mu^2;
          nu   = 2*(2+mu2/v);
    s    = 2*mu*(1+mu2/v);
        out=out+lpdfig2(theta(9),s,nu);
  
         %sigma is inverse Gamma(1,0.5)
   mu=1;
        v=0.5;
        mu2=mu^2;
          nu   = 2*(2+mu2/v);
    s    = 2*mu*(1+mu2/v);
        out=out+lpdfig2(theta(10),s,nu);
        
        
        
        
        
        
        
        
        
