function [ fload,bload ] =...
    getfloadSRtest1( dataS,pmat,FLOAD0,PFLOAD0,...
    rmatin,nfact,NN,...
    rhoin,r0,v0,rho0 )
LI=cols(rhoin);
%draw factor loadings
fload=zeros(NN,nfact);
bload=zeros(NN,1);
rhoout=zeros(NN,LI);
rho0out=rho0;
problemout=0;
for j=1:NN
    yy=dataS(:,j);
    
    xx=[pmat log(rmatin(2:end,j))];
    
    
    %remove serial correlation
    yys=transformrho(yy,rhoin(j,:));
    xxs=transformrho(xx,rhoin(j,:));
    
    %remove heteroscedasticity
    yyh=yys./sqrt(rmatin(2:end,j));
    xxh=xxs./repmat(sqrt(rmatin(2:end,j)),1,cols(xx));
    
    
    yyh=yyh(LI+1:end,:);
    xxh=xxh(LI+1:end,:);
        
        
        fload0=FLOAD0(j,:)';
        pfload0=PFLOAD0;
   
    BB=getreg(yyh,xxh,fload0,pfload0,1);
   
         fload(j,:)=BB(1:nfact)';
         bload(j)=BB(end);
         
         
end
fload(1:nfact,1:nfact)=eye(nfact);

