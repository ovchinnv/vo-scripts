% generate rotations for CHARMM
r=1 ; %radius
np=100000 ; %number of points to represent the spherical shell

% generate points
% NOTE : alternatively, one could generate theta according to sin(phi) and take phi as uniformly distributed
theta=rand(np,1) * 2 * pi ;
phi=acos( - ( rand(np,1) * 2 - 1 ) ) ; % phi is distributed as sin phi (NOTE : can take +/- cos )
zeta=rand(np,1) * 2 * pi ; %a separate  `internal' rotation

% plot to check distribution :
[h,x]=hist(phi);
%plot(x,h/np/(x(2)-x(1)),'k.',x, 0.5*sin(x),'r')
%return
%plot spherical shell :
x=r*sin(phi).*cos(theta);
y=r*sin(phi).*sin(theta);
z=r*cos(phi);

v0=[1 0 0];

rotations = zeros(np,8);

for n=1:np
 v=[x(n), y(n), z(n)]/r;
 vrot=cross(v0,v);
 vn=norm(vrot);
 vrot=vrot/vn;  % rotation vector
 arot=acos(dot(v0,v)); % rotation angle
 rotations(n,:) = [ v0 zeta(n)*180/pi vrot arot*180/pi ]; % two rotations concatenated
end
rout=rotations';rout=rout(:);
save('rotations.dat','-ascii','rout');

