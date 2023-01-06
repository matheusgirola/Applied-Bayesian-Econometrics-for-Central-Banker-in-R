function [out,filter]=hamiltonfilternew(phi1,phi2,sig1,sig2,pr_tr,y,x)

  t=size(y,1);
  n=size(x,2);
  
    %unconditional probabilities
%     
%     ett=[(1-q)/(2-q-p); (1-p)/(2-q-p)];

A = [(eye(2)-pr_tr);ones(1,2)];
           EN=[0;0;1];
           ett11= pinv(A'*A)*A'*EN;


    isig1=1/sig1;
    isig2=1/sig2;
    lik=0;
    filter=zeros(t,2);
    for j=1:t
        em1=y(j)-x(j,:)*phi1; 
        em2=y(j)-x(j,:)*phi2; 
        neta1=(1/sqrt(sig1))*exp(-0.5*(em1*isig1*em1'));%F(Y\S=1)
        neta2=(1/sqrt(sig2))*exp(-0.5*(em2*isig2*em2'));%F(Y\S=2)
        %%%Prediction Step%%%%
        ett10=pr_tr*ett11;
        %%%%Update Step%%%%
        ett11=ett10.*[neta1;neta2]; %joint density F(Y,S)
        fit=sum(ett11);           %Marginal density F(Y)
        ett11=(ett11)/fit;    %conditional density F(S\Y) the weights of the likelihood
        filter(j,1:2)=ett11';      %save filter probability ett  
        if fit>0
            lik=lik+log(fit);
        else
            lik=lik-10;
        end
    end

if isnan(lik) || 1-isreal(lik)
    out=10000000;
else
    out=-lik;
end
