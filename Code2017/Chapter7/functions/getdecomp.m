function [ totalv,totalw,totalc,vole ] = getdecomp( hlastw,hlastc,elast0,beta2w,...
    beta2c,rho,fload,nfact,NC,iamatw,iamatc,Sbigw,Sbigc,L )

T=rows(hlastw)-1;
NN=cols(elast0);
volw=zeros(T,nfact);
volc=zeros(T,nfact*NC);
vole=zeros(T,NN);
totalv=vole;
totalw=vole;
totalc=vole;
fload2=fload.^2;
Li=cols(rho);

for i=1:T
    %vol of world
    sigmaw=iamatw*diag(hlastw(i+1,1).*Sbigw)*iamatw';
    volw(i,:)=volatilityx( beta2w,sigmaw,nfact,L );
    %vol of country
    j=1;
   for jj=1:NC
       sigmac=squeeze(iamatc(jj,:,:))*...
        diag(hlastc(i+1,jj).*Sbigc(jj,:))*squeeze(iamatc(jj,:,:))';
   volc(i,j:j+nfact-1)= volatilityx( beta2c(jj,:),sigmac,nfact,L );
   j=j+nfact;
   end
  %vol of idiosyncratci 
  for jj=1:NN
      sigmae=elast0(i+1,jj);
      vole(i,jj)=volatilityx([rho(jj,:) 0],sigmae,1,Li);
  end
  totalv(i,:)=[volw(i,:) volc(i,:)]*fload2'+vole(i,:);
  totalw(i,:)=[volw(i,:) zeros(1,NC*nfact)]*fload2';
  totalc(i,:)=[zeros(1,nfact) volc(i,:)]*fload2';
end

