function [ fbig,Qbig,mubig,meanvol ] = gettrans(hlast,meanvolin,Qbigin,F0,...
    S0f,g0,vg0,MUF0,SMUF0 )


ytemp=log(hlast(:,1))-meanvolin(1);
xtemp=lag0(ytemp,1);
ytemp=ytemp(2:end,:);
xtemp=xtemp(2:end,:);
[F1,Q1]=getAR(ytemp,xtemp,F0,S0f,Qbigin,g0,vg0);
fbig=F1;
Qbig=Q1;


  tmp=log(hlast(:,1))-lag0(log(hlast(:,1)),1)*fbig;
  ytemp=packr(tmp);
  xtemp=ones(rows(ytemp),1).*(1-fbig);
  ytemp=ytemp./sqrt(Qbig);
  xtemp=xtemp./sqrt(Qbig);
  MUJ=getreg(ytemp,xtemp,MUF0,SMUF0,1);
  mubig=(1-fbig).*MUJ;
  meanvol=MUJ;
end

