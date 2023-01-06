
function[hneww,hnewc,rnew,pmatnew,p0new]= getdensOECDk(b0w,parthw,...
    b0c,parthc,b0r,partr,...
    Sbigw,Sbigc,b0p,p0p,y,MUw,Fw,MUc,Fc,H,MUr,Fr,Qr,...
    nfact,Fx,Mux,iAAw,iAAc,NC)
%draw states

        rtemp=MUr'+b0r.*Fr'+partr.*sqrt(Qr');
        
        rnew=rtemp;
        
        htempw=MUw'+b0w*Fw'+parthw;
        hneww=htempw;

        htempc=MUc'+b0c*Fc'+parthc;

        hnewc=htempc;
       
        %kalman filter
        [pmatnew,p0new]=kfilter(y,H,Fx,Mux,b0p,p0p,...
    rnew,hnewc, hneww,iAAc,iAAw,Sbigc,Sbigw,...
    NC,nfact);