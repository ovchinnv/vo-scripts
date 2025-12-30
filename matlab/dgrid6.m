% custom grid generator -- 4th order polynomial
% copyright 2003 Victor Ovchinnikov
% discrete integration 2005
% symmetry debugging May 2006
% 2025 Freely distributable under MIT license
% NOTE: intentionally produces an extra point beyond the right grid boundary

clear;
close all;
format long;
styles={'k.','r.','g.','b.','m.','c.'};
%
%%%%%%%%%
Nmax=4; % number of resolution specifications ; this usually means there will be N-1 grid segments, because the last spec is the endpoint

Lmax=1; % total grid length
%
% starting grid point for each resolution specification
N(1)=1;
N(2)=N(1)+10;
N(3)=N(2)+100;
N(4)=N(3)+10;

% grid spacing at each specification
V(1)=0.04;
V(2)=0.004;
V(3)=0.004;
V(4)=0.04;
%
% grid slope at each specification
S(1)=-0.03;
S(2)=0.0;
S(3)=0.0;
S(4)=0.03;
%
% distance that each grid specification covers ; for example, for 4 specifications, you need to provide two distances; the third one will be computed using tot length
D(1)=0.3;
D(2)=0.4;
%%%%%%%%%%%
qsmooth=1; %whether to smoothen the final grid ; this is highly recommended
nsmooth=20; % number of smoothing interations
qaddpt=0; % whether to add an extra point to the right of the grid boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (Nmax>2) 
 D(Nmax-1)=Lmax-sum(D(1:Nmax-2));
 else D(Nmax-1)=Lmax;
end;
%
for i=1:Nmax-1;
M(1,:)=[ (N(i)+1)^4, (N(i)+1)^3, (N(i)+1)^2, (N(i)+1), 1]/V(i);
M(2,:)=[ (N(i+1)-2)^4, (N(i+1)-2)^3, (N(i+1)-2)^2, (N(i+1)-2), 1]/V(i+1); % difference between v(N(i+2)-1) & v(N(i+1)-1)

M(3,:)=[N(i)^4,N(i)^3,N(i)^2,N(i),1];
M(4,:)=[(N(i+1)-1)^4,(N(i+1)-1)^3,(N(i+1)-1)^2,(N(i+1)-1),1]; % set spacing at next-to-last point

M(5,:)=[(N(i+1)-N(i))*N(i)^4+4*sum([1:N(i+1)-N(i)-1])*N(i)^3+6*sum([1:N(i+1)-N(i)-1].^2)*N(i)^2+4*sum([1:N(i+1)-N(i)-1].^3)*N(i)+sum([1:N(i+1)-N(i)-1].^4)...
        (N(i+1)-N(i))*N(i)^3+3*sum([1:N(i+1)-N(i)-1])*N(i)^2+3*sum([1:N(i+1)-N(i)-1].^2)*N(i)+sum([1:N(i+1)-N(i)-1].^3)...
        (N(i+1)-N(i))*N(i)^2+2*sum([1:N(i+1)-N(i)-1])*N(i)+sum([1:N(i+1)-N(i)-1].^2)...
        (N(i+1)-N(i))*N(i)+sum([1:N(i+1)-N(i)-1])...
	(N(i+1)-N(i))];


B=[S(i)+1.,1.-S(i+1),V(i),V(i+1),D(i)]';
%=inv(M)*B;
A=M\B;
% now we have the polynomial coefficients in the solution matrix ; use them to compute the grid spacings (g) for each segment
a=A(1);
b=A(2);
c=A(3);
d=A(4);
e=A(5);
 j=N(i);
 g(j)=a*j^4+b*j^3+c*j^2+d*j+e;
 Sum(i)=g(j);
 dxmax(i)=g(j);
 dxmin(i)=g(j);
 for j=N(i)+1:N(i+1)-1;
  g(j)=a*j^4+b*j^3+c*j^2+d*j+e;
  Sum(i)=Sum(i)+g(j);
 end
end

if (qsmooth)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optional filtering to produce a smoother grid %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 x=zeros(1,N(Nmax));
 j=N(1);
 x(j)=0.;
 x(j+1)=g(j);
 for j=N(1)+1:N(Nmax)-1;
  x(j+1)=x(j)+g(j);
 end;
 x=filter1d(x,nsmooth);
 g(1:end)=x(2:end)-x(1:end-1); % recompute the grid spacings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%rescale g to correct length: (e.g. diffusive smoothing above will distort the grid, which will alter its total length)
%sum(Sum)
scale=Lmax/(sum(Sum));
%scale=1.;
Mstr=1.; % maximum stretching
j=N(1);
g(j)=g(j)*scale;
dxmax=g(j);
dxmin=g(j);
lentot=g(j); % running total length
x(j)=0.;
x(j+1)=g(j);
for j=N(1)+1:N(Nmax)-1;
  g(j)=g(j)*scale;
  lentot=lentot+g(j);
  x(j+1)=x(j)+g(j);
  str=g(j)/g(j-1); % compute current stretching coefficient :
  stretch(j)=max(str,1./str); % >1 regardless of whether stretching or compressing
  Mstr=max(Mstr,stretch(j));
  dxmax=max(dxmax,g(j));
  dxmin=min(dxmin,g(j));
end
if (qaddpt)
 x(N(Nmax)+1)=2*x(N(Nmax))-x(N(Nmax)-1); % add extra point
end
stretch(1)=stretch(2);
%
lentot
Mstr
dxmax
dxmin

%fid=fopen('xgrid.dat','w');
%fprintf(fid,'%25.20d\n',x(1:N(Nmax)+1));
%fclose(fid);

figure(1);hold on;
plot(x,'k'); 
title('x vs. np')

figure(2);hold on;
for i=1:Nmax-1 
 plot(linspace(N(i),N(i+1)-1,N(i+1)-N(i)),g(N(i):N(i+1)-1),char(styles(mod(i,6)+1)));hold on;
end
title('dx vs. np');

figure(3);hold on;
plot(x(1:N(Nmax)-1),stretch(1:N(Nmax)-1),'k');
title('stretching factor vs. x');

figure(4);hold on;
for i=1:Nmax-1 
% plot(x(N(i):N(i+1)-1),g(N(i):N(i+1)-1),char(styles(mod(i,6)+1)));hold on;
 plot(x(N(i):N(i+1)-1),g(N(i):N(i+1)-1),char(styles(mod(i,6)+1)));hold on;
end
title('dx vs. x');
xlabel('x');
ylabel('dx'); 
figure(2);
