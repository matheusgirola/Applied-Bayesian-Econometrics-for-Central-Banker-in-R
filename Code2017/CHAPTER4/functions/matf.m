function fm=matf(pm,ns,ps);  


% set initial values for use with iteration @
na = 1;
nb = ns;
nc = ns*ns;
fm = pm;
% iz = 1;
%   do until iz > ps;
  for iz=1:ps
     fz = fm;
     fm = zeros(nc,nc);
%      iw = 1;
%      do until iw > ns;
     for iw=1:ns
         fm(((iw-1)*nb+1):(iw*nb),((iw-1)*na+1):(iw*na))...
            = fz(1:nb,((iw-1)*na+1):iw*na);
     end
%      ib = 2;
%      do until ib > ns;
     for ib=2:ns
        fm(1:nc,((ib-1)*nb+1):ib*nb)  = fm(1:nc,1:nb);
     end
     na = na*ns;
     nb = nb*ns;
     nc = nc*ns;
  end
