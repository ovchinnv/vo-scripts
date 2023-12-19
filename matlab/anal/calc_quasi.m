function [tsquant, tsclass, nu_ps, evc, CX, squant_contrib, sclass_contrib, ev] =calc_quasi(x,y,z,w,temperature);
% perform principal component analysis and quasiharmonic analysis
% note that quasiharmonic frequencies will be incorrect unless mass-weighting is used
%
% some constants :
%
if (~exist('temperature','var'))
 disp(['Temperature was not specified. Will use 300 K ...']);
 temperature=300;
end
%
kb=1.3806503e-23; % J/K
Na=6.0221e23;     %
hbar=1.05457148e-34; % kg m^2 / s   (J . s)
jle=1/4.184/1000;    % one joule in units of kcal
c=2.998e10;        % cm/s in vacuum
kb2=kb*Na*jle;    % kb2 * T = 0.6
kbt=kb2*temperature;
h=hbar*2*pi*Na*jle; %   so that h * freq gives kcal/mol 
tol=1e-10;

%
[natom,nframes]=size(x);
if (~exist('w','var') | isempty(w))
  sw=ones(natom,1);
else
  sw=sqrt(w);
end
sw=repmat(sw,3,1);
%
X=[ x ; y ; z ; ]' ;% concatenate coordinates ; the variables must be in different columns; the observations in different rows
disp(['Computing ',num2str(3*natom),'x',num2str(3*natom),' covariance matrix ...']);
CX=cov(X); % note that this is not the usual procedure ; we usually mass weight the coords first ;
% note also that we are not specifying normalization ; we should use cov(x,1) to nomalize by the # samples N, not by N-1 (default)
% however, for large numbers of samples, the difference will be unimportant
% mean is always subtracted
%
% mass-weight covariance matrix
% diagonalize
disp(['Computing eigenvalues and eigenvectors ...']);
%
[evc,ev]=eig( (sw*sw').*CX); % mass-weight before diagonalization
% scale eigenvectors by inverse square root of mass to obtain Cartesian eigenvectors (modes of vibration)
osw=repmat(1./sw,1,length(sw));
evc=osw.*evc;
%
% compute frequencies
ev=diag(ev); % grab the diagonal entries
% omit below, keeping all evs for information
% NOTE that the RGB modes, if all is well, will be the first 6, but the modes are listed in increasing evalue, or decreasing freq, so low-freq at the end !
%ev=ev(find(ev>tol)); % throw away negative evalues or zeros; this should also remove rigid-body modes (note that they could be at the end ...)
% note that not having enough samples will create additional zero eigenvalues
%
% now have raw eigenvalues from the diagonalization of the covariance matrix (units Ang^2 amu)
% convert frequencies to s^-1 units :
%
fact  = 1./sqrt(jle * 1e-23) / (2*pi);  % 2.0455e13 ! this essentially converts kilocalories to Joules in kT
nu_s  = zeros(size(ev)) ; squant_contrib=nu_s ; sclass_contrib=nu_s ;
okinds = find(ev>tol) ; % conside only positive evs
nu_s(okinds) = sqrt(kbt./ev(okinds)) .* fact ; % per second
nu_ps = nu_s*1e-12 ;% per picosecond
nu_icm= (nu_s/c)';% inverse wavelength per centimeter

% Andricioaei formula
pe=nu_s(okinds)*h/kbt;           % energy of the quasiharmonic modes divided by kBT
%directly
%const = 1./sqrt(kbt   *4.184) * 6.62606 * 6.0221 * 0.1 / 2 / pi;
%pe=     1./sqrt(kbt*ev*4.184) * 6.62606 * 6.0221 * 0.1 / 2 / pi;

squant_contrib(okinds)  =  kb2 * (pe./(exp(pe) - 1) - log (1-exp(-pe))); % this is also the same as the quasiharmonic entropy in numata & knapp
sclass_contrib(okinds)  = -kb2 * (log(pe)-1);

%plot(nu_ps,squant_contrib,'r'); hold on;
%plot(nu_ps,sclass_contrib,'b');

tsquant=temperature*sum(squant_contrib);
tsclass=temperature*sum(sclass_contrib);

return

end
