function out=spectrum(x,omega)

%omega is the desired frequency which can be a vector
%taken from page 167 in Hamilton Time series analysis and uses a Bartlett
%kernel with bandwidth 2*sqrt(T)
T=rows(x);
q=2*sqrt(T);
xbar=mean(x);
j=0;
lamda0=sum((x(j+1:T)-xbar).*(x(j+1:T)-xbar))/T;

out=[];
for jj=1:length(omega)
    f=lamda0;
    f1=0;
    for j=1:q
        xlag=lag0(x,j);
        lamdaj=sum((x(j+1:T)-xbar).*(xlag(j+1:T)-xbar))/T;
        f1=f1+(1-(j/(q+1)))*lamdaj*cos(omega(jj)*j);
    end
    f=(f+2*f1)/(2*pi);
    out=[out;f];
end



















% xbar=mean(x);
% out=[];
% for jj=1:length(J)
% f=0;
% % k=0;
% % lamda0=sum((x(k+1:T)-xbar).*(x(k+1:T)-xbar));
% for k=1:m
%     xlag=lag0(x,k);
%     Rk=sum((x(k+1:T)-xbar).*(xlag(k+1:T)-xbar))/sum((x-xbar).^2);
%     lamdak=1-k/m;
%     omegaj=(J(jj)*(pi))/m;
%     f=f+lamdak*Rk*cos(k*omegaj);
% end
% f=1+2*f;
% 
% out=[out;f/pi];
% end