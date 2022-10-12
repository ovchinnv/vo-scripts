function z=smooth1(x,y,delta)

sx=size(x);
if(sx(1)>1);  xtemp=x'; else ; xtemp=x; end ;

nx=length(x);
z=zeros(1,nx);

for n=1:delta
% correction for external points
  tempx=x(n:delta:nx);%sx=size(tempx);if(sx(1)>1); tempx=tempx';end;
  tempy=y(n:delta:nx);sy=size(tempy);if(sy(1)>1); tempy=tempy';end;
  
  
  if (tempx(1)>x(1))
   tempx=[x(1) tempx(:)'];tempy=[y(1) tempy(:)'];
  end  
  if (tempx(end)<x(end))
   tempx=[tempx(:)' x(end) ];tempy=[tempy(:)' y(end)];
  end  
  
  z=z+interp1(tempx,tempy,xtemp);
end

z=z./delta;
% restore dimension order
if(sx(1)>1);  z=z'; end ;

return;

