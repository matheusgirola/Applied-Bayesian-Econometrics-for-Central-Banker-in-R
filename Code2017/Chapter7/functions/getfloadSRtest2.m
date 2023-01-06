function [ rhoout,problemout,rho0out ]...
    = getfloadSRtest2( dataS,pmat,...
    FLOAD0,PFLOAD0,rmatin,nfact,NN,...
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
    
    %remove heteroscedasticity
    yy=yy./sqrt(rmatin(2:end,j));
    xx=xx./repmat(sqrt(rmatin(2:end,j)),1,cols(xx));
    
   
    BB=FLOAD0(j,:)';     
         %sample rho
     resid=yy-xx*BB;  %no heteroscedasticity but serial correlation
     ye=resid;
     xe=zeros(rows(ye),LI);
     for i=1:LI
         xe(:,i)=lag0(ye,i);
     end
     ye=ye(LI+1:end,:);
     xe=xe(LI+1:end,:);
     [bdraw,problem]=getARx(ye,xe,r0,v0,1);
     if problem
         rhoout(j,:)=rho0(j,:);
     else
          rhoout(j,:)=bdraw';
          rho0out(j,:)=bdraw';
     end
     problemout=problemout+problem;
         
end

