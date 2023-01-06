function [ fload,qout,rhoout,rho0out,problemout,res ] = ...
    getfloadfactggTVP( dataS,pmatw,pmate,pmatc,FLOAD0,PFLOAD0,Q0,Qin,rmatin,nfact,NN,...
    index,NC,FLOAD0e,PFLOAD0e,Q0e,indexe,rhoin,r0,v0,rho0,LI )


T=rows(pmatw);
%draw factor loadings
fload=zeros(NN,T,nfact*3);
rhoout=zeros(NN,LI);
rho0out=rho0;
qout=zeros(NN,nfact*3,nfact*3);
problemout=0;
id=unique(index);
idx=vec(repmat(1:NC,nfact,1));
res=zeros(T,NN);
jx=1;
for ii=1:length(id)
    datax=dataS(:,index==id(ii)); %country i
      NNx=cols(datax);
    pmatx=pmatc(:,idx==id(ii));
    ide=indexe(index==id(ii));
    rhoout1=zeros(NNx,LI);
rho01=rho0(index==id(ii),:);
Qin1=Qin(index==id(ii),:,:);
rho0out1=rho01;
    if ide(1)==1  %EA factor exists
       
    pmat=[pmatw pmate pmatx];
    %emat=eye(cols(pmat));
    fload1=zeros(NNx,T,nfact*3);
    qout1=zeros(NNx,nfact*3,nfact*3);
    NN0=nfact*3;
  FLOAD01=FLOAD0e(index==id(ii),:);
  PFLOAD0x=PFLOAD0e(index==id(ii),:,:);
  Q0x=Q0e(index==id(ii),:,:);
  
  
  
  
  
    else
    pmat=[pmatw pmatx];
   % emat=eye(cols(pmat));
    fload1=zeros(NNx,T,nfact*2);
    qout1=zeros(NNx,nfact*2,nfact*2);
    FLOAD01=FLOAD0(index==id(ii),:,:);
     PFLOAD0x=PFLOAD0(index==id(ii),:,:);
  Q0x=Q0(index==id(ii),:,:);
        NN0=nfact*2;
    end
    
     
   



resx=zeros(T,NNx);
rmatin1=rmatin(:,index==id(ii));
rhoin1=rhoin(index==id(ii),:);



    
for j=1:NNx
    yy=datax(:,j);
    xx=pmat;
    
    %remove serial correlation
    yys=transformrho(yy,rhoin1(j,:));
    xxs=transformrho(xx,rhoin1(j,:));
    %remove heteroscedasticity
    yyh=yys./sqrt(rmatin1(1:end,j));
    xxh=xxs./repmat(sqrt(rmatin1(1:end,j)),1,cols(xx));
    
    
    
    
        
        fload0=FLOAD01(j,:);
        pfload0=squeeze(PFLOAD0x(j,:,:));
       q0=squeeze(Q0x(j,:,:));
       q=squeeze(Qin1(j,1:NN0,1:NN0));
       
       %draw TVP factor loading
       [BB,resf]=carterkohn1(fload0,pfload0,...
            ones(T,1),q,yyh,xxh);

      fload1(j,:,:)=BB;
      %sample q
      resq=diff(BB);
      scale1=resq'*resq+q0;
      T0=cols(q0)+1;
      q1=iwpq(T+T0,invpd(scale1)); %draw from inverse Wishart
      qout1(j,:,:)=q1;
      %sample rho
      resid=yy-sum(xx.*BB,2);  %heteroscedasticity and serial correlation
     ye=resid;
     xe=zeros(rows(ye),LI);
    
     for i=1:LI
         xe(:,i)=lag0(ye,i);
     end
     yeh=ye./sqrt(rmatin1(1:end,j));
    xeh=xe./repmat(sqrt(rmatin1(1:end,j)),1,cols(xe));
     
     
     yeh=yeh(LI+1:end,:);
     xeh=xeh(LI+1:end,:);
     [bdraw,problem]=getARx(yeh,xeh,r0,v0,1);
     if problem
         rhoout1(j,:)=rho01(j,:);
     else
          rhoout1(j,:)=bdraw';
          rho0out1(j,:)=bdraw';
     end   
      problemout=problemout+problem;
    %%%%%%%%%%%remove serial correlation to save residuals for SVOL step
    yyss=transformrho(yy,rhoout1(j,:));
    xxss=transformrho(xx,rhoout1(j,:)); 
         
    resx(:,j)=yyss-sum(xxss.*BB,2); %residuals have no serial correlation by heteroscedastic (epsilon)
end
 
if ide(1)==1
      fload1x=fload1;
 else
    fload1x=zeros(NNx,T,nfact*3);
    fload1x(:,:,1:nfact)=fload1(:,:,1:nfact);
    fload1x(:,:,(nfact*2)+1:nfact*3)=fload1(:,:,nfact+1:end);
 end

fload(jx:jx+NNx-1,:,:)=fload1x;
qout(jx:jx+NNx-1,1:NN0,1:NN0)=qout1;
res(:,jx:jx+NNx-1)=resx;
rhoout(jx:jx+NNx-1,:)=rhoout1;
rho0out(jx:jx+NNx-1,:)=rho0out1;
jx=jx+NNx;
end