function out=pchy(x)
lnx=log(x);
lnxlag=lag0(lnx,4);
out=(lnx-lnxlag)*100; %annual growth