function [ fload,bload ] = getfloadSL( dataS,pmat,FLOAD0,PFLOAD0,rmatin,nfact,NN )
%draw factor loadings
fload=zeros(NN,nfact);
bload=zeros(NN,size(rmatin,3));

for j=1:NN
    yy=dataS(:,j);
%     size(pmat)
%     size(squeeze(log(rmatin(j,2:end,:))))
    xx=[pmat squeeze(log(rmatin(j,2:end,:)))];
    
    yy=yy./sqrt(rmatin(j,2:end,1))';
    xx=xx./repmat(sqrt(rmatin(j,2:end,1))',1,cols(xx));
        
        
        fload0=FLOAD0(j,:)';
        pfload0=PFLOAD0;
   
    BB=getreg(yy,xx,fload0,pfload0,1);
   
         fload(j,:)=BB(1:nfact)';
         bload(j,:)=BB(nfact+1:end);
     
end
fload(1:nfact,1:nfact)=eye(nfact);

