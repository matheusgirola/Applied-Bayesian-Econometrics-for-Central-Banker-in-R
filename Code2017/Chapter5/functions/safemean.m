function pqm=safemean(pfload)
sfactor=max(pfload);
temp=exp(pfload-sfactor);
pqm=nanmean(temp);
pqm=log(pqm)+sfactor;