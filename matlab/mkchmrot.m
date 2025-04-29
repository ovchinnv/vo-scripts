% generate rotations for CHARMM ; based on ~/idea/rotation ; use with urs-rot.inp
% here, I am determining rotations as the effect of another rotation on aprescribed vector, e.g. [1 0 0]
if exist('OCTAVE_VERSION') ; graphics_toolkit('gnuplot') ; end
r=1 ; %radius

if ~exist('qcheck') ; qcheck = 0 ; end
if ~exist('qrand')  ; qrand = 1  ; end
if ~exist('qinternal')  ; qinternal = 1  ; end

if (qrand)
 np=200 ; %number of points to represent the spherical shell
% generate points
% NOTE : alternatively, one could generate theta according to sin(phi) and take phi as uniformly distributed
 theta=rand(np,1) * 2 * pi ;
 phi=acos( - ( rand(np,1) * 2 - 1 ) ) ; % phi is distributed as sin phi (NOTE : can take +/- cos )
% plot to check phi distribution :
%[h,x]=hist(phi);
%plot(x,h/np/(x(2)-x(1)),'k.',x, 0.5*sin(x),'r')
 if (qinternal)
  nzeta=np;
  zeta=rand(nzeta,1) * 2 * pi ; %a separate  `internal' rotation
 end
 x=r*sin(phi).*cos(theta);
 y=r*sin(phi).*sin(theta);
 z=r*cos(phi);
 ntot=np;
else
 nphi=10 ;
 nth=10 ;
 nzeta=2 ;
 ival=0; % initial value for angle distributions ; setting this to 1 removes the many 0 rotations ; but then the distribution is slightly skewed from uniformity
% the skewness improves as you increase the number of points
% generate points
% exclude one endpoint :
 th = ( 2*pi/nth * [ ival : nth-1 ] ) ; % since theta is uniformly distributed, any point has equal weight, so generate a uniform mesh
% phi=acos( - ( linspace(0,1,nphi) * 2 - 1 ) ) ; % phi is distributed as sin phi ; to get the transformation function, integrate & invert, get -acos 
 phi=acos( - ( linspace(0,1,nphi+1) * 2 - 1 ) ) ; phi=phi(1:nphi); % not sure whether it is better to exclude the last point ... seems to make no diff ...
% to exclude the 0th phu angle :
 if (ival==1); phi=phi(2:end) ; nphi=nphi-1; end
 zeta = ( 2*pi/nzeta * [ ival : nzeta-1 ] ) ; % for "internal" rotation
% you can think of linspace as a kind of sort(rand(0,1)), so the sample should be correctly distributed
 [TH,PHI]=ndgrid(th,phi);
 x=r*sin(PHI(:)).*cos(TH(:)); % distribution  uniform for large points
 y=r*sin(PHI(:)).*sin(TH(:));
 z=r*cos(PHI(:)); % exactly uniform (observed, and by construction above)
 np=numel(x) ;
 ntot=np*max(1,(nzeta-ival)*qinternal);
end

if (qcheck)
%plot spherical shell :
 figure(1) ; clf
 ms=15;
 scatter3(x(:),y(:),z(:),ms,'k');
end

v0=[1 0 0];

rotations = zeros(ntot,4*(qinternal+1));

for nz = 1 : ntot/np
for n=1:np
 v=[x(n), y(n), z(n)]/r; % to make unit if not already (r=1)
 vrot=cross(v0,v);
 vn=norm(vrot);
 vrot=vrot/vn;  % rotation vector
 arot=acos(dot(v0,v)); % rotation angle
 if (qinternal)
  rotations(np*(nz-1)+n,:) = [ v0 zeta(nz)*180/pi vrot arot*180/pi ]; % two rotations concatenated
 else
  rotations(n,:) = [ vrot arot*180/pi ];
 end
end %n
end %nz
rout=rotations';rout=rout(:);
if (qrand)
 save(['rotations-rand-',num2str(ntot),'.dat'],'-ascii','rout');
else
 save(['rotations-grid-',num2str(ntot),'.dat'],'-ascii','rout');
% save(['rotations-grid-',num2str(ntot),'.dat'],'-ascii','rotations');
end
