function [ outL,outh,fevdL,fevdh ] = getirfsign(L,LH,N,horizon1,Fmat,Qmat,varcoef,A0  )
mse=0;
for i=1:N
    shock=zeros(1,N);
    shock(i)=1;
      tmp= getlevelirf( L,LH,N,horizon1,shock,Fmat,Qmat,varcoef,A0 );
       tmpv= getvolirfs( L,LH,N,horizon1,shock,Fmat,Qmat,varcoef );
       outL(i,:,:)=tmp;
       outh(i,:,:)=tmpv;
       mse=mse+(cumsum(tmp.^2)+cumsum(tmpv.^2));
       
      
end

for i=1:N
    shock=zeros(1,N);
    shock(i)=1;
      tmp= getlevelirf( L,LH,N,horizon1,shock,Fmat,Qmat,varcoef,A0 );
       tmpv= getvolirfs( L,LH,N,horizon1,shock,Fmat,Qmat,varcoef );
       fevdL(i,:,:)=cumsum(tmp.^2)./mse;
       fevdh(i,:,:)=cumsum(tmpv.^2)./mse;

end

