function [ A ] = formA0( Theta )
[R,C]=size(Theta);
if C==1
  A=diag(Theta(1:3));
    
    A(2,1)=Theta(4);
    A(3,2)=Theta(5);
else
    A=[diag(Theta);Theta(2,1);Theta(3,2)];
end


end

