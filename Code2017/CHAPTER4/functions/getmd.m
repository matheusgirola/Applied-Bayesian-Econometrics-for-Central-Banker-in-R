function [m,v]=getmd(alpha)
salpha=sum(alpha);
m=alpha./salpha;
v=(alpha.*(salpha-alpha))./((salpha.^2)*(salpha+1));