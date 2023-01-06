function [y,x] = create_dummiesSSVOL(lamda,tau,delta,epsilon,p,mu,sigma,n,ph,epsilonH,sigmaH,TR)
% Creates matrices of dummy observations [...];
%lamda tightness parameter
%tau  prior on sum of coefficients
%delta prior mean for VAR coefficients
% epsilon tigtness of the prior around constant
% mu sample mean of the data
% sigma AR residual variances for the data




% Initialise output (necessary for final concatenation to work when tau=0):
x = [];
y = [];
yd1 = [];
yd2 = [];
xd1 = [];
xd2 = [];
nn=1;
%% Get dummy matrices in equation (5) of Banbura et al. 2007:

	yd1=[diag(sigma.*delta)./lamda;
         zeros(n*(p-1),n)];
     for i=1:ph*nn
         yd1=[yd1;zeros(1,n)];
     end
     
         if TR==0
        yd1=[yd1; zeros(1,n)];
         else
                yd1=[yd1; zeros(2,n)];
         end     
     
	jp=diag(1:p);
    jph=diag(1:ph);
%     kron(jph,1./sqrt(epsilonH))
	xd1=[kron(jp,diag(sigma)./lamda)];
    xdh1=kron(jph,diag(sigmaH)./epsilonH);
    if TR==0
    tmp=blkdiag(xd1,xdh1,epsilon);
    else
    tmp=blkdiag(xd1,xdh1,epsilon,epsilon);
    end

%     for i=1:ph*nn
%          xd1=[[xd1 zeros(rows(xd1),1)]; [zeros(1,cols(xd1)) sigmaH./sqrt(epsilonH)]];
%      end
%    
%          xd1=[[xd1 zeros(rows(xd1),1)]; [zeros(1,cols(xd1)) epsilon]];
   

%% Get additional dummy matrices - see equation (9) of Banbura et al. 2007:
if tau>0
    
	yd2=diag(delta.*mu)./tau;
    if TR==0
	xd2=[kron(ones(1,p),yd2) zeros(n,nn*ph) zeros(n,1)];
    else
     	xd2=[kron(ones(1,p),yd2) zeros(n,nn*ph) zeros(n,2)];
    end
   
    
end
     
%% 
y=[yd1;yd2];
% xd1
% size(xd1)
% size(xd2)
% size(kron(ones(1,p),yd2))
size(tmp)
size(xd2)
x=[tmp;xd2];
 
         
 
 