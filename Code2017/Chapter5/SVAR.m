clear;
addpath('functions','sims_Optimization');
bx1=[ 0.7 -0.1  0.1 0;  -0.1 0.7 0.1 0; -0.1 0.1 0.7 0 ];
A0x=[0.1 0  0;
    -0.1 0.1 0;
    0   0.5  0.1];
sigmax=inv(A0x)*inv(A0x)';



T=500;
N=3;
dataout=zeros(T,N);

for i=2:T  
    dataout(i,:)=[dataout(i-1,:) 1]*bx1'+randn(1,3)*chol(sigmax);
end
irftrue1=irfsim(bx1',N,1,inv(A0x)',[1 0 0],21);
irftrue2=irfsim(bx1',N,1,inv(A0x)',[0 1 0],21);
irftrue3=irfsim(bx1',N,1,inv(A0x)',[0 0 1],20);

data=dataout;
L=1;

Y=data;
N=cols(Y);
ncrit=(N*L+1);

%take lags
X=[];
for j=1:L
X=[X lag0(data,j) ];
end
X=[X ones(rows(X),1)];

B=X\Y;
iXX=inv(X'*X);
E=Y-X*B;
SB=E'*E;
T=rows(X);
K=cols(X);
theta0=ones(5,1).*0.1;
out=getML(theta0,SB,T,K);
options = optimset('Disp','iter','Diagnostics','on','LargeScale','off',...
    'MaxFunEvals',100000,'MaxIter',500,'TolFun',1e-05,'TolX',1e-05);
% 
% 
% % %simplex
 [Theta1,fval] = fminsearch(@getML, theta0,options,SB,T,K);
 % %sims
 [fval,Theta2,gh,hess,itct,fcount,retcodeh] = csminwel('getML',Theta1,eye(length(Theta1))*.5,[],1e-15,1000,SB,T,K);
 
 A0ML=formA0(Theta2);
 %check all diagonal elements are positive and switch sign if not
 indx=find(diag(A0ML)<0);
 A0ML(:,indx)=A0ML(:,indx)*-1;
 
 %%%%%%%%%%%%%%MCMC Algorithm%%%%%%%%%%%%%%%%%%%%%%%%%%%
 REPS=30000;
 BURN=20000;
 thetaold=formA0(A0ML);
 scale=0.6;
 P=chol(hess)*scale;
 naccept=0;
 jj=1;
 for j=1:REPS
     %draw from candidate
     thetanew=thetaold+(randn(1,rows(thetaold))*P)';
     
     %acceptance probability
     liknew=-getML( thetanew,SB,T,K );
     likold=-getML( thetaold,SB,T,K );
     
      accept=liknew-likold;
    
    if accept>log(rand)
        thetaold=thetanew;
        naccept=naccept+1;
    end
    arate=naccept/j;
    
    %Normalise A0
     A0=formA0(thetanew);
      indx=find(diag(A0)<0);
 A0(:,indx)=A0(:,indx)*-1;
% indx=find(diag(A0\A0ML)<0); %alternative normalisation used by Waggoner
% Zha
% A0(:,indx)=A0(:,indx)*-1;
    %draw coefficients conditional on A0
   
    sigma=inv(A0)*inv(A0)';
    V=kron(sigma,iXX);
    beta=vec(B)+(randn(1,N*(N*L+1))*chol(V))';
    
    disp(sprintf(' Replication %s of %s.', ... 
             num2str([j arate] ), num2str(REPS)) );
    
    if j>BURN;
        %compute IRFs
 irf1=irfsim(reshape(beta,N*L+1,N),N,L,inv(A0)',[1 0 0],21);
irf2=irfsim(reshape(beta,N*L+1,N),N,L,inv(A0)',[0 1 0],21);
irf3=irfsim(reshape(beta,N*L+1,N),N,L,inv(A0)',[0 0 1],21);

out1(jj,:,:)=irf1;
out2(jj,:,:)=irf2;
out3(jj,:,:)=irf3;
jj=jj+1;
    end
 end
 
 %plot
 tmp1=prctile(out1,[50 16 84]);
 tmp2=prctile(out2,[50 16 84]);
 tmp3=prctile(out3,[50 16 84]);
 figure(1)
 for j=1:3
     subplot(3,3,j);
     plot(tmp1(:,:,j)','r');
     hold on
     plot(irftrue1(:,j),'k');
 end
  jj=4;  
 for j=1:3
     subplot(3,3,jj);
     plot(tmp2(:,:,j)','r');
     hold on
     plot(irftrue2(:,j),'k');
     jj=jj+1;
 end   
    
   jj=7;  
 for j=1:3
     subplot(3,3,jj);
     plot(tmp3(:,:,j)','r');
     hold on
     plot(irftrue3(:,j),'k');
     jj=jj+1;
 end     
 
