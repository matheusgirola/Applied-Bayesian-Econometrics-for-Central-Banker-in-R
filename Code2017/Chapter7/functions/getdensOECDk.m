
function[hneww,hnewc,rnew,pmatnew,p0new,floadnew]= getdensOECDk(b0w,parthw,b0e,parthe,...
    b0c,parthc,b0r,partr,...
    Sbigw,Sbige,Sbigc,b0p,p0p,y,MUw,Fw,MUe,Fe,MUc,Fc,MUr,Fr,Qr,...
    nfact,Fx,Mux,iAAw,iAAe,iAAc,NC,b0f,parthf,NN,NSf,rho,index,id,L)

%draw states
        floadnew=zeros(NN,NSf);
        for ii=1:NN
            ftemp=(b0f(ii,:))+parthf(NN,:);
            floadnew(ii,:)=ftemp;
        end
     
        rtemp=MUr'+b0r.*Fr'+partr.*sqrt(Qr');
        
        rnew=rtemp;
        
        htempw=MUw'+b0w*Fw'+parthw;
        hneww=htempw;
        
          htempe=MUe'+b0e*Fe'+parthe;
        hnewe=htempe;

        htempc=MUc'+b0c*Fc'+parthc;

        hnewc=htempc;
       
        %kalman filter
        [pmatnew,p0new]=kfilter(y,floadnew,Fx,Mux,b0p,p0p,...
    rnew,hnewc, hneww,hnewe,iAAc,iAAw,iAAe,Sbigc,Sbigw,Sbige,...
    NC,nfact,rho,index,id,L);