function out=woodbury(Ainv,C,Binv)
	
	%inv(A+C*B*C')
  
%  	out=Ainv-Ainv*C*invpd(Binv+C'*Ainv*C)*C'*Ainv;
      AinvC=Ainv*C;
      Ct=C';
	out=Ainv-AinvC*invpd(Binv+Ct*AinvC)*Ct*Ainv;