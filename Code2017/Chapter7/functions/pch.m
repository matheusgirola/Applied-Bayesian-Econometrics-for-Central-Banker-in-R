function out=pch(x)
lnx=log(x);
lnxlag=lag0(lnx,1);
out=(lnx-lnxlag)*400; %quarterly annualised growth