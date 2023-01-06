clear
addpath('./functions');
dfolder='./data/';
sfolder='./results/';
file=1;
sfile=strcat(dfolder,'dataxx0',num2str(file));
%Load data and transform%%%%%%%
load(sfile);
dataS=standardise(dataS); %standardise data
%%estimation options%%%%%%
T0=20;   %training sample
L=2;     %lag for transition equation
Lx=1;    %lag for idiosyncratic component transition eq
REPS=20000; %Reps
BURN=10000; %burn-in
SKIP=10;  %every SKIP draw is kept after burn-in 
maxdraws=100; %max trys to find stable coefficients
CHECK=1;
Sindex=BURN+1:SKIP:REPS;
fsize=length(Sindex);
id=unique(index); %index of countries
NC=length(id); %number of countries
NN=cols(dataS); %number of series
idc=vec(repmat(1:NC,1,1)); %index of countries

%%%%Starting Values and Priors
%initial estimate of the factors
pmatw=extract(dataS,1); %PC estimator of world factor
dataSS=dataS-pmatw*(pmatw\dataS);
pmatc=zeros(rows(pmatw),NC); %PC for countries

for i=1:NC
    dataC=dataSS(:,index==id(i));
    tmp=extract(dataC,1);
pmatc(:,i)=tmp;

end

res=zeros(rows(pmatw),NN); %idiosyncratic
FLOAD0=zeros(NN,2);  %prior mean for Factor loadings
for j=1:NN
      
        yy=dataS(:,j);
        xx=[pmatw pmatc(:,idc==index(j))];
        BB=xx\yy;
        res(:,j)=yy-xx*BB;
        FLOAD0(j,:)=BB';
end
VFLOAD0=eye(2).*10; %prior variance

%priors for TVP parameters
scale=3.5e-04;
%world factor
[y0w,x0w]=preparex(pmatw(1:T0,:),L,1);
[b00w,s00w,p00w]=getols(y0w,x0w); %OLS AR regression on pre-sample
Q0w=scale*p00w*T0; %OLS covariance times T0 times scaling (scale matrix for IW prior Qw~IW(Q0w, T0)
Qw=Q0w; %starting value
%country factors
b00c=cell(NC,1);
s00c=cell(NC,1);
p00c=cell(NC,1);
Q0c=cell(NC,1);
Qc=zeros(L+1,L+1,NC);

for j=1:NC
[y0c,x0c]=preparex(pmatc(1:T0,j),L,1);
[b00c{j},s00c{j},p00c{j}]=getols(y0c,x0c);
Q0c{j}=scale*p00c{j}*T0;  %Qc~IW(Q0c{j},T0) for j=1,2,...NC
Qc(:,:,j)=scale*p00c{j}*T0; %starting value for Qc

end
%idiosyncratic
b00e=cell(NN,1);
s00e=cell(NN,1);
p00e=cell(NN,1);
Q0e=cell(NN,1);
Qe=zeros(Lx,Lx,NN);
for j=1:NN
[y0e,x0e]=preparex(res(1:T0,j),Lx,0);
[b00e{j},s00e{j},p00e{j}]=getols(y0e,x0e);
Q0e{j}=scale*p00e{j}*T0; %Qe~IW(Q0e{j},T0) for j=1,2,...NN
Qe(:,:,j)=scale*p00e{j}*T0; %Starting values
end

%remove training sample
dataS=dataS(T0+1:end,:);
pmatw=pmatw(T0+1:end,:);
pmatc=pmatc(T0+1:end,:);
res=res(T0+1:end,:);
T=rows(dataS);


%priors and starting values for stochastic volatilties as residual^2+small
%number
%world
[y0w,x0w]=preparex(pmatw,L,1);
[~,~,~,epsw]=getols(y0w,x0w); %regression of factor on lags
hlastw=epsw.^2+0.0001; %residual^2+small number
hlastw=[hlastw(1:L+1,:);hlastw];
%country
hlastc=zeros(T+1,NC);

for j=1:NC
[y0c,x0c]=preparex(pmatc(:,j),L,1);
[~,~,~,epsc]=getols(y0c,x0c);%regression of factor on lags
hlastcc=epsc.^2+0.0001;%residual^2+small number
hlastcc=[hlastcc(1:L+1,:);hlastcc];
hlastc(:,j)=hlastcc;
end

%idiosyncratic
hlaste=zeros(T+1,NN);
for j=1:NN
[y0e,x0e]=preparex(res(:,j),Lx,0);
[~,~,~,epse]=getols(y0e,x0e);%regression of factor on lags
hlastee=epse.^2+0.0001;%residual^2+small number
hlastee=[hlastee(1:Lx+1,:);hlastee];
hlaste(:,j)=hlastee;
end

SS0=10;    %variance of initial condition of SVOL
g0=0.1^2;  %prior scale parameter for inverse gamma prior for g
Tg0=1;     %prior degrees of freedom
gw=g0;  %starting values
gc=ones(NC,1).*g0;  %starting values
ge=ones(NN,1).*g0;   %starting values

beta0w=repmat(b00w',T,1);
beta0c=zeros(T,L+1,NC);
for j=1:NC
    beta0c(:,:,j)=repmat(b00c{j}',T,1);
end
beta0e=zeros(T,Lx,NN);
for j=1:NN
    beta0e(:,:,j)=repmat(b00e{j}',T,1);
end

%initial conditions for the factors
pmat00=[pmatw(L,:) pmatc(L,:)];
for j=1:L-1
    pmat00=[pmat00 [pmatw(L-j,:) pmatc(L-j,:)]];
end
vmat00=eye(cols(pmat00))*1;

save priors

decompsave=zeros(fsize,T,NN);
pmatsave=zeros(fsize,T,NC+1);
hsave=zeros(fsize,T+1,NC+1);

jgibbs=1;
igibbs=1;

while jgibbs<=fsize
%%%%%%%%%%%%%%%%%Gibbs Step 1: Draw TVP Parameters%%%%%%%%%%%%%%%%%
%1 A: World Factor 
[yw,xw]=preparex([pmatw(1:L,:);pmatw],L,1);
[beta2w,errorw,rootsw,problemw]=...
    carterkohnAR(yw,xw,Qw,hlastw,b00w',p00w,L,CHECK,maxdraws,1);
if problemw
    beta2w=beta0w;
else
    beta0w=beta2w;
end
%draw Qw
resbeta=diff(beta2w);
scaleQ=resbeta'*resbeta+Q0w;
Qw=iwpq(T+T0,invpd(scaleQ));
%1 B: Country Factors
beta2c=zeros(T,L+1,NC);
errorc=zeros(T,NC);
problemC=zeros(NC,1);
for j=1:NC
    [yc,xc]=preparex([pmatc(1:L,j);pmatc(:,j)],L,1);
    [beta2c(:,:,j),errorc(:,j),rootsc,problemc]=...
    carterkohnAR(yc,xc,Qc(:,:,j),hlastc(:,j),b00c{j}',p00c{j},L,CHECK,maxdraws,1);
if problemc
    beta2c(:,:,j)=beta0c(:,:,j);
else
    beta0c(:,:,j)=beta2c(:,:,j);
end
%draw Qc
resbeta=diff(beta2c(:,:,j));
scaleQ=resbeta'*resbeta+Q0c{j};
Qc(:,:,j)=iwpq(T+T0,invpd(scaleQ));
problemC(j)=problemc;
end
%1 c: Idiosyncratic Factors
beta2e=zeros(T,Lx,NN);
errore=zeros(T,NN);
problemE=zeros(NN,1);
parfor j=1:NN
    [ye,xe]=preparex([res(1:Lx,j);res(:,j)],Lx,0);
    [beta2e(:,:,j),errore(:,j),rootse,probleme]=...
    carterkohnAR(ye,xe,Qe(:,:,j),hlaste(:,j),b00e{j}',p00e{j},Lx,CHECK,maxdraws,0);
if probleme
    beta2e(:,:,j)=beta0e(:,:,j);
else
    beta0e(:,:,j)=beta2e(:,:,j);
end
%draw Qe
resbeta=diff(beta2e(:,:,j));
scaleQ=resbeta'*resbeta+Q0e{j};
Qe(:,:,j)=iwpq(T+T0,invpd(scaleQ));
problemE(j)=probleme;
end

%%%%%%%%%%%%%%%%%Gibbs Step 2: Draw SVOL%%%%%%%%%%%%%%%%%    
%%2A: World factor SVOL
hneww=getsvol(hlastw,gw,log(s00w),SS0,errorw);
hlastw=hneww;
gerrors=diff(log(hlastw));
gw=IG(Tg0,g0,gerrors);
%%2B: Country factor SVOL
hnewc=zeros(T+1,NC);
for j=1:NC
    hnewc(:,j)=getsvol(hlastc(:,j),gc(j),log(s00c{j}),SS0,errorc(:,j));
    hlastc(:,j)=hnewc(:,j);
    gerrors=diff(log(hlastc(:,j)));
    gc(j)=IG(Tg0,g0,gerrors);
end
 %% 2C: Idiosyncratic factor SVOL  
  hnewe=zeros(T+1,NN);
parfor j=1:NN
    hnewe(:,j)=getsvol(hlaste(:,j),ge(j),log(s00e{j}),SS0,errore(:,j));
    hlaste(:,j)=hnewe(:,j);
    gerrors=diff(log(hlaste(:,j)));
    ge(j)=IG(Tg0,g0,gerrors);
end  

%% Step 3: Sample factor loadings

fload=zeros(NN,2);
res=zeros(T,NN);
jjj=1;
for jj=1:NC
 tmpdata=dataS(:,idc(jj)==index);     %data for country i
 tmpbeta2e=beta2e(:,:,idc(jj)==index);
 tmphlaste=hlaste(:,idc(jj)==index);
 tmpfload0=FLOAD0(idc(jj)==index,:);
 fload1=zeros(cols(tmpdata),2);
 res1=zeros(T,cols(tmpdata));
 fload1(1:2,:)=eye(2); %identification
 for j=1:2
     yy=tmpdata(:,j);
        xx=[pmatw pmatc(:,jj)];
     res1(:,j)=yy-xx*fload1(j,:)';
 end
      for j=3:cols(tmpdata)
        yy=tmpdata(:,j);
        xx=[pmatw pmatc(:,jj)];
        %remove serial correlation
        yys=remSC(yy,tmpbeta2e(:,:,j));
        xxs=remSC(xx,tmpbeta2e(:,:,j));
        %remove heteroscedasticity
        yyss=yys./sqrt(tmphlaste(2:end,j));
        xxss=xxs./repmat(sqrt(tmphlaste(2:end,j)),1,cols(xxs));
        %take care of missing values
        yyss=yyss(Lx+1:end,:);
        xxss=xxss(Lx+1:end,:);
        %draw from conditional posterior
        FL=getreg(yyss,xxss,tmpfload0(j,:)',VFLOAD0,1);
        fload1(j,:)=FL';
        %save residuals
        res1(:,j)=yy-xx*FL; %residuals are serially correlated and heteroscedastic
      end
    fload(jjj:jjj+cols(tmpdata)-1,:)=fload1; %save factor loadings for each country
    res(:,jjj:jjj+cols(tmpdata)-1)=res1;
    jjj=jjj+cols(tmpdata);    
end

%% Step 4: Carter Kohn Algorithm to sample the factors
dataF=zeros(T,NN);
for j=1:NN
    dataF(:,j)=remSC(dataS(:,j),beta2e(:,:,j));
end
dataF(1:Lx,:)=repmat(dataF(Lx+1,:),Lx,1);
%Carter and Kohn algorithm to draw the factor
ns=cols(pmat00);
beta_tt=zeros(T,ns);          %will hold the filtered state variable
ptt=zeros(T,ns,ns);    % will hold its variance
beta11=pmat00;
p11=vmat00;

for i=1:T
    
    %build matrices of state space as they are time-varying
    
    %observation equation
    H1=zeros(NN,NC+1);
    H2=H1;
    %world factor loadings
    H1(:,1)=fload(:,1);
    H2(:,1)=fload(:,1).*-squeeze(beta2e(i,:,:));
    %country factor loadings
    jj=2;
    jjj=1;
    for j=1:NC
        floadc=fload(index==idc(j),2); %country factor loadings
        tmpbeta2e=squeeze(beta2e(i,:,idc(j)==index)); %AR coefficient at time t of idiosyncratic shock
     H1(jjj:jjj+rows(floadc)-1,jj)=floadc;   
     H2(jjj:jjj+rows(floadc)-1,jj)=floadc.*-tmpbeta2e;  
     jj=jj+1;
     jjj=jjj+rows(floadc);
    end
    H=zeros(NN,(NC+1)*2);
    H(:,1:NC+1)=H1;
    H(:,NC+2:end)=H2;
    
    R=diag(hlaste(i+1,:));
    
    %transition equation
    
    Q=zeros(ns,ns);
    Q(1:NC+1,1:NC+1)=diag([hlastw(i+1) hlastc(i+1,:)]);
    F1=diag([beta2w(i,1) ;squeeze(beta2c(i,1,:))]); %AR 1 coefficients
    F2=diag([beta2w(i,2) ;squeeze(beta2c(i,2,:))]); % AR 2 coefficients
    F=[[F1 F2];eye(ns-(NC+1),ns)];
    MU=[beta2w(i,L+1)  squeeze(beta2c(i,L+1,:))' zeros(1,ns-(NC+1))];
    
    %Prediction
x=H;
beta10=MU+beta11*F';
p10=F*p11*F'+Q;
yhat=(x*(beta10)')';                                               
eta= dataF(i,:)-yhat;
feta=(x*p10*x')+R;
ifeta=invpd(feta);
%updating
K=(p10*x')*ifeta;
beta11=(beta10'+K*eta')';
p11=p10-K*(x*p10);
ptt(i,:,:)=p11;
beta_tt(i,:)=beta11;

end
% Backward recursion to calculate the mean and variance of the distribution of the state
%vector
beta2 = zeros(T,ns);   %this will hold the draw of the state variable
jv1=1:NC+1; %index of state variables to extract
jv=jv1;
wa=randn(T,ns);

i=T;  %period t
p00=squeeze(ptt(i,jv1,jv1)); 
beta2(i,:)=beta_tt(i,:);
beta2(i,jv1)=beta_tt(i:i,jv1)+(wa(i:i,jv1)*cholx(p00));   %draw for beta in period t from N(beta_tt,ptt)
%periods t-1..to .1
for i=T-1:-1:1

 %build matrices of transition equation

    Q=zeros(ns,ns);
    Q(1:NC+1,1:NC+1)=diag([hlastw(i+2) hlastc(i+2,:)]);
    F1=diag([beta2w(i+1,1) ;squeeze(beta2c(i+1,1,:))]); %AR 1 coefficients
    F2=diag([beta2w(i+1,2) ;squeeze(beta2c(i+1,2,:))]); % AR 2 coefficients
    F=[[F1 F2];eye(ns-(NC+1),ns)];
    MU=[beta2w(i+1,L+1)  squeeze(beta2c(i+1,L+1,:))' zeros(1,ns-(NC+1))];
    
f=F(jv,:);
q=Q(jv,jv);
mu=MU(jv);
pt=squeeze(ptt(i,:,:));
ifptfq=invpd(f*pt*f'+q);
bm=beta_tt(i:i,:)+(pt*f'*ifptfq*(beta2(i+1:i+1,jv)-mu-beta_tt(i,:)*f')')';  
pm=pt-pt*f'*ifptfq*f*pt;  
beta2(i,:)=bm;
beta2(i:i,jv1)=bm(jv1)+(wa(i:i,jv1)*cholx(pm(jv1,jv1)));  
end


pmat=beta2(:,jv1);   %update the factors
pmatw=pmat(:,1);
pmatc=pmat(:,2:NC+1);

if igibbs>BURN
    if sum(Sindex==igibbs)>0
%% calculate variance decomposition

% Volatility of the world factor
VOLW=zeros(T,1);
for i=1:T
    VOLW(i,:)=volatility(beta2w(i,:)',hlastw(i+1,:),1,L);
end

%Volatility of the country factor
VOLC=zeros(T,NC);
for j=1:NC
    for i=1:T
        VOLC(i,j)=volatility(beta2c(i,:,j)',hlastc(i+1,j),1,L);
    end
end

%volatility of the idiosyncratic factor
VOLE=zeros(T,NN);
for j=1:NN
    for i=1:T
        VOLE(i,j)=volatility([beta2e(i,:,j)';0],hlaste(i+1,j),1,Lx);
    end
end

floadsquared=H1.^2;

 TOTALVOL=[VOLW VOLC]*floadsquared'+VOLE;      
  TOTALVOLW=[VOLW zeros(rows(VOLC),cols(VOLC))]*floadsquared';     
  DECOMP=TOTALVOLW./TOTALVOL;
  
  decompsave(jgibbs,:,:)=DECOMP;
  pmatsave(jgibbs,:,:)=[pmatw pmatc];
  hsave(jgibbs,:,:)=[hlastw hlastc];
  jgibbs=jgibbs+1;
    end
end

disp(strcat('REPS=',num2str([igibbs ])));
igibbs=igibbs+1;
end

%% Plot results

tmp=prctile((pmatsave),[50 16 84]);

tt=1:T;
figure(1)
for j=1:5
    subplot(2,3,j)
    [hh,h1,h2]=plotyy(tt,(tmp(:,:,j)'),tt,pmat00x(21:end,j));
    set(hh(1),'xlim',[min(tt) max(tt)])
    set(hh(2),'xlim',[min(tt) max(tt)])
end









tmp=prctile(log(hsave),[50 16 84]);

tt=1:T+1;
figure(2)
for j=1:5
    subplot(2,3,j)
    [hh,h1,h2]=plotyy(tt,(tmp(:,:,j)'),tt,hlast0x(20:end,j));
    set(hh(1),'xlim',[min(tt) max(tt)])
    set(hh(2),'xlim',[min(tt) max(tt)])
end

    
    tmp=prctile(decompsave,[50 ]);

figure(3)
for j=1:40
    subplot(10,4,j)
    plot(tmp(:,:,j)'*100,'r')
    hold on
    plot(vdecompwx(21:end,j)*100,'k')
    title(strcat('\pi_{',num2str(j),',t}'),'Interpreter','tex')
    ylim([0 100])
end

figure(4)
for j=1:40
    subplot(10,4,j)
    plot(tmp(:,:,j+40)'*100,'r')
    hold on
    plot(vdecompwx(21:end,j+40)*100,'k')
    title(strcat('\pi_{',num2str(j+40),',t}'),'Interpreter','tex')
ylim([0 100])
end