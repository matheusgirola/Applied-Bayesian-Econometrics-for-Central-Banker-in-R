function [ ir1,ir2,ir3,ir4] = get_irfTVAR4NOTRENDSTAR(N, HORZ, T, L,LH, Y0,hlast0, beta1, Sbig1,iamat1,beta2,tar,tvar,delay,muvol,fvol,qvol,reps,scale,nml,EX,gam)



% Draw N(0,1) innovations for variance  equation:
Cq=sqrt(qvol);
bb1=reshape(beta1,N*L+EX,N);
bb2=reshape(beta2,N*L+EX,N);
% Note we only need HORZ*N innovations, but adding an extra L draws makes 
% the indexing in the loop below a bit cleaner.
yy1=0;
yy2=0;
yy3=0;yy4=0;yy=0;
for ii=1:reps
%compute forecast
yhat=zeros(HORZ+L,N);
yhat(1:L,:)=Y0;
hhat=zeros(HORZ+L,1);
hhat(L-LH+1:L,:)=log(hlast0);

yhat1=yhat;
yhat2=yhat;
yhat3=yhat;
yhat4=yhat;
ystar=zeros(HORZ+L,1);

ystar4=ystar;
ystar1=ystar;
ystar2=ystar;
ystar3=ystar;
for fi=L+1:HORZ+L
    hhat(fi,:)=muvol+hhat(fi-1,:)*fvol+randn(1,1)*Cq;
    
    xhat=[];
    xhat1=[];
    xhat2=[];
    xhat3=[];
    xhat4=[];
    for ji=1:L
        xhat=[xhat yhat(fi-ji,:)];
        xhat1=[xhat1 yhat1(fi-ji,:)];
        xhat2=[xhat2 yhat2(fi-ji,:)];
        xhat3=[xhat3 yhat3(fi-ji,:)];
        xhat4=[xhat4 yhat4(fi-ji,:)];
    end
    
    for ji=0:LH
        xhat=[xhat hhat(fi-ji,:)];
        xhat1=[xhat1 hhat(fi-ji,:)];
        xhat2=[xhat2 hhat(fi-ji,:)];
        xhat3=[xhat3 hhat(fi-ji,:)];
        xhat4=[xhat4 hhat(fi-ji,:)];
    end
    
    xhat=[xhat 1 ];
    xhat4=[xhat4 1];
    xhat1=[xhat1 1 ];
    xhat2=[xhat2 1 ];
    xhat3=[xhat3 1 ];
  
    ystar(fi,:)=yhat(fi-delay,tvar);
    ystar4(fi,:)=yhat4(fi-delay,tvar);
    ystar1(fi,:)=yhat1(fi-delay,tvar);
    ystar2(fi,:)=yhat2(fi-delay,tvar);
    ystar3(fi,:)=yhat3(fi-delay,tvar);
    
    %build covariance matrix
    sigma1=iamat1*diag(Sbig1.*exp(hhat(fi,1)))*iamat1';
    csigma1=cholx(sigma1);
    d1=diag(csigma1);
    csx1=csigma1./repmat(d1,1,N);
    
    
    %simulate data
    if fi==L+1
    uu1=zeros(1,N);
           
    e2=1./(1+exp(-gam.*(ystar(fi,:)-tar)));
    
    yhat(fi,:) = (xhat*bb1)+((xhat*bb2)*e2)+uu1*csigma1;
    
    
    
     %
     %y4
    uu1=zeros(1,N);uu1(4)=scale;
     e2=1./(1+exp(-gam.*(ystar4(fi,:)-tar)));
    if nml==0
    yhat4(fi,:) = (xhat4*bb1)+((xhat4*bb2)*e2)+uu1*csigma1;
    else
        yhat4(fi,:) = (xhat4*bb1)+((xhat4*bb2)*e2)+uu1*csx1;
    end
    
    %y1
    uu1=zeros(1,N);uu1(1)=scale;
     e2=1./(1+exp(-gam.*(ystar1(fi,:)-tar)));
    if nml==0
    yhat1(fi,:) = (xhat1*bb1)+((xhat1*bb2)*e2)+uu1*csigma1;
    else
    yhat1(fi,:) = (xhat1*bb1)+((xhat1*bb2)*e2)+uu1*csx1;
    end
    
    
  %y2
    uu1=zeros(1,N);uu1(2)=scale;
     e2=1./(1+exp(-gam.*(ystar2(fi,:)-tar)));
    if nml==0
    yhat2(fi,:) = (xhat2*bb1)+((xhat2*bb2)*e2)+uu1*csigma1;
    else
    yhat2(fi,:) = (xhat2*bb1)+((xhat2*bb2)*e2)+uu1*csx1;
    end
    
    
    %y3
    uu1=zeros(1,N);uu1(3)=scale;
     e2=1./(1+exp(-gam.*(ystar3(fi,:)-tar)));
    if nml==0
    yhat3(fi,:) = (xhat3*bb1)+((xhat3*bb2)*e2)+uu1*csigma1;
    else
    yhat3(fi,:) = (xhat3*bb1)+((xhat3*bb2)*e2)+uu1*csx1;
    end
    
    
    
    else
    uu1=randn(1,N);      
    e2=1./(1+exp(-gam.*(ystar(fi,:)-tar)));
    yhat(fi,:) = (xhat*bb1)+((xhat*bb2)*e2)+uu1*csigma1;
    
     %
    uu1=randn(1,N);      
    e2=1./(1+exp(-gam.*(ystar4(fi,:)-tar)));
    yhat4(fi,:) = (xhat4*bb1)+((xhat4*bb2)*e2)+uu1*csigma1;
    
    %
        uu1=randn(1,N);      
    e2=1./(1+exp(-gam.*(ystar1(fi,:)-tar)));
    yhat1(fi,:) = (xhat1*bb1)+((xhat1*bb2)*e2)+uu1*csigma1;
    
    %
        uu1=randn(1,N);      
    e2=1./(1+exp(-gam.*(ystar2(fi,:)-tar)));
    yhat2(fi,:) = (xhat2*bb1)+((xhat2*bb2)*e2)+uu1*csigma1;
    %
       uu1=randn(1,N);      
    e2=1./(1+exp(-gam.*(ystar3(fi,:)-tar)));
    yhat3(fi,:) = (xhat3*bb1)+((xhat3*bb2)*e2)+uu1*csigma1;
    
    end
   
end

yy=yy+yhat;
   yy1=yy1+yhat1;
   yy2=yy2+yhat2;
   yy3=yy3+yhat3;
   yy4=yy4+yhat4;

end

yy=yy/reps;
yy1=yy1/reps;
yy2=yy2/reps;
yy3=yy3/reps;
yy4=yy4/reps;


ir4=yy4-yy;
ir1=yy1-yy;
ir2=yy2-yy;
ir3=yy3-yy;

 ir4=ir4(L+1:end,:);
 ir1=ir1(L+1:end,:);
 ir2=ir2(L+1:end,:);
 ir3=ir3(L+1:end,:);
