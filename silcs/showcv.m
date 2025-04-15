% plot cvs from history file to make sure their range is consistent with the restraints
%
if exist('OCTAVE_VERSION') ; graphics_toolkit('gnuplot') ; end

ngroups=8;
numcv=nchoosek(ngroups,2)+ngroups;

d=load('hist.dat');
d=reshape(d,numcv,[])
% ind1 = 1, ngroups+1, 2*ngroups, 3*ngroups-2, 4*ngroups-5 ...
ind1=[1] ; i=1 ; while ngroups>numel(ind1) ; i=i+ngroups-numel(ind1)+1 ; ind1=[ind1 i]; end
ind2=setdiff(1:numcv,ind1)
d1=d(ind1,:); % spherical restraints
d2=d(ind2,:); % distance restraints

figure(1)
plot(d1') ; % this shows that the restraint values are consistent -- what else could be wrong ?

figure(2) ;
plot(d2') ;