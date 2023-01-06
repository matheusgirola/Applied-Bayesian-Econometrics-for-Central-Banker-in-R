function [ out] = getML( Theta,SB,T,K )
%check bounds
    %calculate lik
  A=formA0(Theta);
    
    lastterm=trace(A*SB*A');
    lik=(T-K)*log(det(A))-(1/2)*lastterm;
    
    if isinf(lik) || ~isreal(lik) || isnan(lik)
        out=inf;
    else
        out=-lik;
    end





