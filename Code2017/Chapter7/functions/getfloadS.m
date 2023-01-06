function [ fload,bload ] = getfloadS( dataS,pmat,FLOAD0,PFLOAD0,rmatin,nfact,NN )
%draw factor loadings
fload=zeros(NN,nfact);
bload=zeros(NN,1);

for j=1:NN
    yy=dataS(:,j);
    
    xx=[pmat log(rmatin(2:end,j))];
    
    yy=yy./sqrt(rmatin(2:end,j));
    xx=xx./repmat(sqrt(rmatin(2:end,j)),1,cols(xx));
        
        
        fload0=FLOAD0(j,:)';
        pfload0=PFLOAD0;
   
    BB=getreg(yy,xx,fload0,pfload0,1);
   
         fload(j,:)=BB(1:nfact)';
         bload(j)=BB(end);
     
end
fload(1:nfact,1:nfact)=eye(nfact);

