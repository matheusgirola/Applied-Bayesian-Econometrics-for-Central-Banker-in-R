function [hnew,naccept]=svol1C(mu0,iP00,hleadx,FF,QQ,iQQ,mumat)
NS=cols(mu0);
XX=log(mu0);
hlead=hleadx;
BINV=iP00+(FF'*iQQ*FF);
b=(iP00*XX')'+(log(hlead)-mumat')*iQQ*FF;


VV=invpd(BINV);
MM=VV*b';
htrial=(exp(MM+(randn(1,NS)*chol(VV))'))';
    hnew=htrial;
    naccept=1;
