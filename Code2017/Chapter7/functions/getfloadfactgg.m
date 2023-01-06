function [ fload,problemout,res ] = ...
    getfloadfactgg( dataS,pmatw,pmate,pmatc,FLOAD0,PFLOAD0,PFLOAD01,rmatin,nfact,NN,...
    index,NC,FLOAD0e,PFLOAD0e,PFLOAD01e,indexe )
maxtrys=100; %hard coded

T=rows(pmatw);
%draw factor loadings
fload=zeros(NN,nfact*3);

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
    if ide(1)==1  %EA factor exists
       
    pmat=[pmatw pmate pmatx];
    emat=eye(cols(pmat));
    fload1=zeros(NNx,nfact*3);
    NN0=nfact*3;
  FLOAD01=FLOAD0e(index==id(ii),:);
  PFLOAD01x=PFLOAD01e;
  PFLOAD0x=PFLOAD0e;
    else
    pmat=[pmatw pmatx];
    emat=eye(cols(pmat));
    fload1=zeros(NNx,nfact*2);
    FLOAD01=FLOAD0(index==id(ii),:);
     PFLOAD01x=PFLOAD01;
       PFLOAD0x=PFLOAD0;
        NN0=nfact*2;
    end
    
     
   



resx=zeros(T,NNx);
rmatin1=rmatin(:,index==id(ii));

for j=1:NN0
    yy=datax(:,j);
    xx=pmat(:,1:j);
    
   
    
    %remove heteroscedasticity
    yyh=yy./sqrt(rmatin1(1:end,j));
    xxh=xx./repmat(sqrt(rmatin1(1:end,j)),1,cols(xx));
    

        
        
        fload0=FLOAD01(j,1:cols(xxh))';
        fload0=fload0*sign(fload0(j)); %make sure diagonal positive 
        pfload0=PFLOAD01x(1:cols(xxh),1:cols(xxh));
        problemc=0;
    check=-1;
    trys=1;
    while check<0 && trys<=maxtrys
       
    BB=getreg(yyh,xxh,fload0,pfload0,1);
%     [fload0 BB]
    if BB(j)>0
        check=10;
    else
       
        trys=trys+1;
    end
    end
    if check<0
        BB=fload0;
        problemc=1;
    end
    BBtemp=emat(j,:)';
    BBtemp(1:cols(xxh),1)=BB;
   
    fload1(j,:)=BBtemp';
    
    
    resx(:,j)=yy-xx*BB; %with heteroscedasticity
 problemout=problemout+problemc;        
end

    
for j=(NN0)+1:NNx
    yy=datax(:,j);
    xx=pmat;
    
    
    %remove heteroscedasticity
    yyh=yy./sqrt(rmatin1(1:end,j));
    xxh=xx./repmat(sqrt(rmatin1(1:end,j)),1,cols(xx));
  
        
        fload0=FLOAD01(j,:)';
        pfload0=PFLOAD0x;
   
    BB=getreg(yyh,xxh,fload0,pfload0,1);
   
         fload1(j,:)=BB';
         
         
         
    resx(:,j)=yy-xx*BB;
end
 if ide(1)==1
      fload1x=fload1;
 else
    
     fload1x=[fload1(:,1:nfact) zeros(NNx,nfact) fload1(:,nfact+1:end)]; %zeros for OECD  factor loadings when NA
 end


fload(jx:jx+NNx-1,:)=fload1x;

res(:,jx:jx+NNx-1)=resx;
jx=jx+NNx;
end