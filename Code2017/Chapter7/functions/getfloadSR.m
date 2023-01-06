function [ fload,bload,rhoout,problemout,rho0out ] = getfloadSR( dataS,pmat,FLOAD0,PFLOAD0,rmatin,nfact,NN,...
    rhoin,r0,v0,rho0 )
LI=cols(rhoin);
%draw factor loadings
fload=zeros(NN,nfact);
bload=zeros(NN,1);
rhoout=zeros(NN,LI);
rho0out=rho0;
problemout=0;
emat=eye(nfact);
for j=1:nfact
    yy=dataS(:,j)-pmat*emat(j,:)';
    
    xx=log(rmatin(2:end,j));
    
    %remove serial correlation
    yys=transformrho(yy,rhoin(j,:));
    xxs=transformrho(xx,rhoin(j,:));
    
    %remove heteroscedasticity
    yyh=yys./sqrt(rmatin(2:end,j));
    xxh=xxs./repmat(sqrt(rmatin(2:end,j)),1,cols(xx));
    
    
    yyh=yyh(LI+1:end,:);
    xxh=xxh(LI+1:end,:);
        
        
        fload0=FLOAD0(j,nfact+1:end)';
        pfload0=PFLOAD0(nfact+1:end,nfact+1:end);
   
    BB1=getreg(yyh,xxh,fload0,pfload0,1);
    BB=[emat(j,:)';BB1];
   
         fload(j,:)=BB(1:nfact)';
         bload(j)=BB(end);
         
         %sample rho
     resid=yy-xx*BB1;  %heteroscedasticity and serial correlation
     ye=resid;
     xe=zeros(rows(ye),LI);
    
     for i=1:LI
         xe(:,i)=lag0(ye,i);
     end
     yeh=ye./sqrt(rmatin(2:end,j));
    xeh=xe./repmat(sqrt(rmatin(2:end,j)),1,cols(xe));
     
     
     yeh=yeh(LI+1:end,:);
     xeh=xeh(LI+1:end,:);
     [bdraw,problem]=getARx(yeh,xeh,r0,v0,1);
     if problem
         rhoout(j,:)=rho0(j,:);
     else
          rhoout(j,:)=bdraw';
          rho0out(j,:)=bdraw';
     end
     problemout=problemout+problem;
         
end
    
for j=nfact+1:NN
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
         
         %sample rho
     resid=yy-xx*BB;  %heteroscedasticity and serial correlation
     ye=resid;
     xe=zeros(rows(ye),LI);
    
     for i=1:LI
         xe(:,i)=lag0(ye,i);
     end
     yeh=ye./sqrt(rmatin(2:end,j));
    xeh=xe./repmat(sqrt(rmatin(2:end,j)),1,cols(xe));
     
     
     yeh=yeh(LI+1:end,:);
     xeh=xeh(LI+1:end,:);
     [bdraw,problem]=getARx(yeh,xeh,r0,v0,1);
     if problem
         rhoout(j,:)=rho0(j,:);
     else
          rhoout(j,:)=bdraw';
          rho0out(j,:)=bdraw';
     end
     problemout=problemout+problem;
         
end

