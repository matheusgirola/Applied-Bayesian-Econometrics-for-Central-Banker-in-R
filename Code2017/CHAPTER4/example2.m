clear;
addpath('functions');
%generate artificial data from MS Model
T=500;
B1=0.2;
B2=0.9;
C1=1;
C2=-1;
S1=3;
S2=1;
P=[0.95 0.05;0.05 0.95];
strue=zeros(T,2);
strue(1,1)=1; %initial state
strue=simS(strue,P); %generate state variable
e=randn(T,1);
Y=zeros(T,1);
X=zeros(T,1);
for i=2:T;
    X(i,:)=Y(i-1,:);
    if strue(i,1)==1
    Y(i)=[X(i,:) 1]*[B1 C1]'+e(i)*sqrt(S1);
    else
    Y(i)=[X(i,:) 1]*[B2 C2]'+e(i)*sqrt(S2);
    end
end



 %%%%%%%%%%%%%%%%Run Hamilton Filter%%%%%%%%%%%%%%%%
   %unconditional probabilities

A = [(eye(2)-P);ones(1,2)];
           EN=[0;0;1];
           ett11= pinv(A'*A)*A'*EN;
    iS1=1/S1;
    iS2=1/S2;
    lik=0;
    filter=zeros(T,2);
    for j=1:T
        em1=Y(j)-[X(j,:) 1]*[B1 C1]'; 
        em2=Y(j)-[X(j,:) 1]*[B2 C2]'; 
        neta1=(1/sqrt(S1))*exp(-0.5*(em1*iS1*em1'));%F(Y\S=1)
        neta2=(1/sqrt(S2))*exp(-0.5*(em2*iS2*em2'));%F(Y\S=2)
        %%%Prediction Step%%%%
        ett10=P*ett11;
        %%%%Update Step%%%%
        ett11=ett10.*[neta1;neta2]; %joint density F(Y,S)
        fit=sum(ett11);           %Marginal density F(Y)
        ett11=(ett11)/fit;    %conditional density F(S\Y) the weights of the likelihood
        filter(j,1:2)=ett11';      %save filter probability ett  
        lik=lik+log(fit);      %save log likelihood
        
    end
  %%%%%%%%Backward Recursion to draw S%%%%%%  
    
  S=zeros(T,1);
   %time T
   p1=filter(T,1);
   p2=filter(T,2);
   p=p1/(p1+p2);
   u=rand(1,1);
   S(T,1)=(u>=p);
  
   for t=T-1:-1:1
   if S(t+1)==0
p00=P(1,1)*filter(t,1);
p01=P(1,2)*filter(t,2);
elseif S(t+1)==1
p00=P(2,1)*filter(t,1);
p01=P(2,2)*filter(t,2);
   end
  u=rand(1,1);
  p=p00/(p00+p01);
  if u<p
      S(t)=0;
  else
      S(t)=1;
  end
   end



