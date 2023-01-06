clear

%generate artificial data
T=100;
sigma=1;
b1=4;
b2=2;
e=randn(T,1)*sqrt(sigma);
X=rand(T,1);
Y=b1.*(X.^b2)+e;



%step 1 set starting values

Gammaold=[0;0;0.1];

%step 2 set SIGMA matrix via OLS estimation
yols=Y;
xols=[ones(T,1) X];
bols=inv(xols'*xols)*(xols'*yols);
eols=yols-xols*bols;
sols=((eols'*eols)/T);
vols=sols*inv(xols'*xols);

P=eye(3);                    %this is the variance of the metropolis hastings random walk based partly on OLS estimates
P(1,1)=(vols(1,1));
P(2,2)=(vols(2,2));
P(3,3)=0.1;
K=0.2;
P=K*P;

REPS=15000;
true=repmat([b1 b2 sigma],REPS,1);   

% step 3 metropolis Hastings algorithm 
out=[];
naccept=0;
for j=1:REPS

    %step 3a draw new Gamma
    Gammanew=Gammaold+(randn(1,3)*chol(P))';
    %step 3b evaluate posterior at new draw
    B1=Gammanew(1);B2=Gammanew(2);sigma2=Gammanew(3);
    if sigma2<0
        posteriorNEW=-1000000;
    else
    B=[B1;B2];
    lik=-(T/2)*log(2*(pi))-(T/2)*sigma2-0.5*(((Y-B1.*(X.^B2))'*(Y-B1.*(X.^B2)))/sigma2); %likelihood function
    posteriorNEW=lik; %posterior at the new draw
    end
    %step 3c evaluate posterior at new draw
    B1=Gammaold(1);B2=Gammaold(2);sigma2=Gammaold(3);
    B=[B1;B2];
    lik=-(T/2)*log(2*(22/7))-(T/2)*sigma2-0.5*(((Y-B1.*(X.^B2))'*(Y-B1.*(X.^B2)))/sigma2); %likelihood function
    posteriorOLD=lik; %posterior at the old draw
    
    %step 3d compute acceptance probability
    accept=min([exp(posteriorNEW-posteriorOLD);1]);   %min(accept,1)
    
    u=rand(1,1);  %random number from the uniform dist
    
    if u<accept
        Gammaold=Gammanew;  %accept draw
        naccept=naccept+1;  %count number of acceptances
    end
    out=[out;Gammaold'];
    
end
 
plot([out true])

xlabel('Metropolis Hastings Draws');
legend('B_1','B_2','\sigma^{2}','True B_1','True B_2','True \sigma^{2}');
