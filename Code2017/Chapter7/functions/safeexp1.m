function [ temp,sfactor] = safeexp1( pfload )

sfactor=max(pfload);
pfload+sfactor
temp=exp(pfload-sfactor);


end