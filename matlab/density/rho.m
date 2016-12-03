% compute density
setup;

dz=0.1 ; % histogram resolution in Ang

savename='rho_dmpc.mat' ;

istep=[];
sx=[];
sy=[];
sz=[];

% get box sizes first
% note that I somehow have truncated xst data for last (12) run ...
for fname = fnames
 xst = load([char(fname),'.xstfile']);
 istep=[istep ; xst(:,1) ];
 sx=[sx ; xst(:,2) ];
 sy=[sy ; xst(:,6) ];
 sz=[sz ; xst(:,10) ];
end
% number of frames
nframes=length(istep) ;
% determine number of histogram bins
zmax=max(sz/2) ; % assume symmetry
nbins=floor(zmax/dz)+1;
zbins=([1:nbins]-0.5)*dz ;
ed=zeros(nbins,nframes) ; % electron density matrix

% go over trajectories to reach nframes
iframe=0 ;

for fname = fnames
 dcdfile = [char(fname),'.dcd'];
 disp(['==> Processing file ',dcdfile]);
 h=read_dcdheader(dcdfile) ;
% loop over frames
 for k=1:h.NSET
  [x,y,z]=read_dcdstep(h) ;
  iframe=iframe+1;
  disp(['==> Processing frame #', num2str(iframe)]);
  if (iframe>nframes) ; break ; end ;% break if number of frames exceeded
%
% translate to center of membrane
% NOTE: if membrane is shifted far, it might be wrapped, so be careful (e.g. look at it in VMD first)
% repeat twice :
  for l=1:2
   zcen = sum( z(lipid)'.*mass(lipid) ) / sum(mass(lipid)) ;
   z=z-zcen ; % make zero correspond to center of membrane
% wrap z coordinates : 
   zz=sz(iframe) ; % cell size
%   z=mod(z+zz/2,zz) ; % wrap, but shift into positive region for histogramming
   z=mod(z+zz/2,zz)-zz/2 ; % wrap
  end
% now histogram :
%  ibin=floor(z/dz) + 1 ; % bin index 
  ibin=floor(abs(z)/dz) + 1 ; % bin index 
% see no other way but to loop over atoms
  for j=1:natom
   ed(ibin(j),iframe) = ed(ibin(j),iframe) + anumber(j) ;
  end
   ed(:,iframe) = ed(:,iframe) / sx(iframe) / sy(iframe) / dz / 2 ;
 end
end 

iave=floor(0.5*nframes) ; % average only over latter half of data
%
eda=mean(ed(:,iave:end),2); % average density
save(savename,'ed','zbins','istep');

% smooth average : 
addpath '~/scripts/matlab/' ;
edas=smooth2(zbins,eda,8);

plot(zbins,edas) ; 
xlim([0 36]);

% location of maximum

[mx,imx] = max(edas) ;
zbins(imx)

