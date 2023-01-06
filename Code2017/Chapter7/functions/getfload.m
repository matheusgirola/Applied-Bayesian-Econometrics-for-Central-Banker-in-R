function [ fload,rmatnew ] = getfload( dataS,pmat,FLOAD0,PFLOAD0,rmatin,nfact,NN,TF0,VF0 )
%draw factor loadings
fload=zeros(NN,nfact);
rmatnew=zeros(NN,1);
for j=1:NN
    yy=dataS(:,j);
    
        xx=pmat;
        fload0=FLOAD0(j,:)';
        pfload0=PFLOAD0;
   
    BB=getreg(yy,xx,fload0,pfload0,rmatin(j));
    resf=yy-xx*BB;
     rmati= IG(TF0,VF0,resf);
     rmatnew(j)=rmati;
     
         fload(j,:)=BB';
     
end
fload(1:nfact,1:nfact)=eye(nfact);

