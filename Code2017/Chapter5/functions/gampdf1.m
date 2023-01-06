function Y=gampdf1(v,delta,h)
X=h;
A=v/2;
B=2/delta;

Y = log(gampdf(X,A,B));