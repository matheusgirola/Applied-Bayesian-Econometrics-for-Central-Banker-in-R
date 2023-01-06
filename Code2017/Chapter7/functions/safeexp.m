function [ temp,sfactor] = safeexp( pfload )
sfactor=max(pfload);
temp=exp(pfload-sfactor);


end