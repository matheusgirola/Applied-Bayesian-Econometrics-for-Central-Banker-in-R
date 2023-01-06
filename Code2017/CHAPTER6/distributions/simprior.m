function out=simprior(ndraws)
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





out=[];
%sigma has a gamma prior with amean of 1 and variance 1
%compute shape and scale parameters
        mu=1;
        v=1;
        b = v/mu;
        a = mu/b;
        outi=gamrnd(a,b,ndraws,1);
        out=[out outi];
        
 %delta has to be Gamma prior

        
          mu=1.5;
        v=1;
        b = v/mu;
        a = mu/b;
        outi=gamrnd(a,b,ndraws,1);
        out=[out outi];
         
  %alpha has gamma prior with mean 3 and variance 1
  
        mu=3;
        v=1;
        b = v/mu;
        a = mu/b;
        outi=gamrnd(a,b,ndraws,1);;%lpdfgam(theta(3),a,b);
        out=[out outi];
  
  
  %omega is Gamma 
      mu=1.5;
        v=1;
        b = v/mu;
        a = mu/b;
        outi=gamrnd(a,b,ndraws,1);
        out=[out outi];
   
   
   %rho is beta with mean 0.5 variance 0.2
    mu=0.5;
     v=0.2;
     a = (1-mu)*mu^2/v^2 - mu;
        b = a*(1/mu - 1);
        out =[out betarnd(a,b,ndraws,1)];% lpdfbeta(theta(5),a,b);
        
        
        %rho is beta with mean 0.5 variance 0.2
     mu=0.5;
     v=0.2;
     a = (1-mu)*mu^2/v^2 - mu;
        b = a*(1/mu - 1);
        out =[out betarnd(a,b,ndraws,1)];
        
        
        %rho is beta with mean 0.5 variance 0.2
     mu=0.5;
     v=0.2;
     a = (1-mu)*mu^2/v^2 - mu;
        b = a*(1/mu - 1);
        out =[out betarnd(a,b,ndraws,1)];
        
        
        %sigma is inverse Gamma(1,1)
         mu=1;
        v=0.5;
        mu2=mu^2;
          nu   = 2*(2+mu2/v);
    s    = 2*mu*(1+mu2/v);
      a=(nu/2);
      b=(2/s);
        out=[out 1./gamrnd(a,b,ndraws,1)]; %see Equation A.18 pg 292 Bauwens et.al
        
  
   %sigma is inverse Gamma(1,1)
        mu=1;
        v=0.5;
        mu2=mu^2;
          nu   = 2*(2+mu2/v);
    s    = 2*mu*(1+mu2/v);
      a=(nu/2);
      b=(2/s);
        out=[out 1./gamrnd(a,b,ndraws,1)]; %see Equation A.18 pg 292 Bauwens et.al
  
         %sigma is inverse Gamma(1,1)
       mu=1;
        v=0.5;
        mu2=mu^2;
          nu   = 2*(2+mu2/v);
    s    = 2*mu*(1+mu2/v);
      a=(nu/2);
      b=(2/s);
        out=[out 1./gamrnd(a,b,ndraws,1)]; %see Equation A.18 pg 292 Bauwens et.al
  
        
        
        
        
        
        
        
        
        
        
        
        
        
        

