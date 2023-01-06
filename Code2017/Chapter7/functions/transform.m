function out=transform(x,id)
 out=[];
 for i=1:size(x,2)
     if id(i)==1;
         temp=diff(log(x(:,i)))*1200;
     elseif id(i)==2
         temp=x(2:end,i)*12;
     elseif id(i)==0
         temp=x(1:end,i);
     elseif id(i)==5
         temp=log(x(2:end,i))*100;
     elseif id(i)==6;
         temp=hpfilter(log(x(1:end,i)),14400)*100;
 
     elseif id(i)==7
         temp=diff(log(x(:,i)))*400;
         elseif id(i)==8
         temp=x(2:end,i)*4;
           elseif id(i)==9;
         temp=hpfilter(x(:,i),14400);
     
     end
     out=[out temp];
 end
 end

