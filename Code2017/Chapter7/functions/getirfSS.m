function [ outL,outh,fevdL,fevdh ] = getirfSS(L,LH,N,horizon1,Fmat,Qmat,varcoef,A0  )
mse=0;
for i=1:N
    shock=zeros(1,N);
    shock(i)=1;
      tmp= getlevelirfSS( L,LH,N,horizon1,shock,Fmat,Qmat,varcoef,A0 );
       outL(i,:,:)=tmp;
      mse=mse+(cumsum(tmp.^2));
       
      
end

       tmpv= getvolirfSS( L,LH,N,horizon1,1,Fmat,Qmat,varcoef );
 outh=tmpv;
       mse=mse+(cumsum(tmpv.^2));

for i=1:N
    shock=zeros(1,N);
    shock(i)=1;
      tmp= getlevelirfSS( L,LH,N,horizon1,shock,Fmat,Qmat,varcoef,A0 );
       fevdL(i,:,:)=cumsum(tmp.^2)./mse;
      

end
       tmpv= getvolirfSS( L,LH,N,horizon1,1,Fmat,Qmat,varcoef );
 fevdh=cumsum(tmpv.^2)./mse;
