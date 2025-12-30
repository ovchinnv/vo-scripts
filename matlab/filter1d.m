% 1-dimensional tophat filtering
function [output]=filter1d(x,n)

 l=length(x);
 for i=1:n
   x(2:l-1)=0.25*(x(1:l-2)+2*x(2:l-1)+x(3:l));
 end;
 output=x;
 return;
