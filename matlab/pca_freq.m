% plot quasiharmonic frequencies 

pf=load('ppsf.dat');
pf=pf(:,2:2:end);

[n,m]=size(pf); nfreq=m*n;
pf=reshape(pf',nfreq,1);
pf=pf(7:end);

rf=load('rigf.dat');
rf=rf(:,2:2:end);

[n,m]=size(rf); nfreq=m*n;
rf=reshape(rf',nfreq,1);
rf=rf(7:end);

% compute entropy
%
kb=1.3806503e-23; % J/K
Na=6.0221e23;     %
hbar=1.05457148e-34; % kg m^2 / s   (J . s)
jle=1/4.184/1000;      % in kcal
T=300;            % K
c=2.998e10;        % cm/s in vacuum
kb2=kb*Na*jle;    % kb2 * T = 0.6
kbt=kb2*T;
h=hbar*2*pi*Na*jle; %   so that h * pf2 gives kcal/mol 

%%%%%%% starting at 2000 gives the number cvlose to what CHARMM spits out
ibeg=1;
pf2=c.*pf(ibeg:end);   %/2/pi;     % 1/second
pe=pf2*h/kbt;           % energy of the quasiharmonic modes divided by kBT
ps= kb2* sum (pe./(exp(pe) - 1) - log (1-exp(-pe)))  %entropy of PPS
%
rf2=c.*rf(ibeg:end);   %/2/pi;     % 1/second
re=rf2*h/kbt;           % energy of the quasiharmonic modes divided by kBT
rs= kb2* sum (re./(exp(re) - 1) - log (1-exp(-re)))  %entropy of RIG

(ps-rs)*T

