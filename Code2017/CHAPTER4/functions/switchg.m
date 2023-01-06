function out=switchg(s,g)
     

     n = rows(s);
     m = rows(g);
     swt = zeros(m,m);     

     t = 2;
     while t<=n;

     st1 = s(t-1);
     st = s(t);
     swt(st1,st) = swt(st1,st) + 1;
     t = t+ 1;
     end

     out=swt;