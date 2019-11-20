% MATLAB routines for trajectory analysis
% coarse-grain by putting average residue coordinates into CA
%
addpath('~/scripts/matlab/anal');
name='traf-rank';

pdbfile=['/home/taly/vmd/prepare/traf/rank/traf-rank.pdb']; % in lieu of a "structure" file

struc ;
% define some selections
select ;

% dcd trajectory names :
basename=['wsh7-ndt-rf/',name];
flag='-now';
ext='dcd';
%
inds=1 ;
dcdnames={};
for i = inds
 dcdnames=[dcdnames { [basename,flag,'.',ext] } ];
end
%
% read dcds :
traj;
%
% initialize CG coordinates
xpdbcg=zeros(natom,1);
ypdbcg=zeros(natom,1);
zpdbcg=zeros(natom,1);
%
xcg=zeros(natom,nall);
ycg=zeros(natom,nall);
zcg=zeros(natom,nall);
%
% coarse grain system to one bead per residue
% search for all CA atoms; identify them by resid and insertion code
%
disp(['Coarse-graining system to one (CA) bead per residue...']);
cainds=find(typeCA);
ncg=length(cainds);

for caind = cainds'
% resinds = find( resid==resid(caind) & insertion==insertion(caind) & ismember(segid,segid(caind))); % make sure to match the segid also
 resinds = find( resid==resid(caind) & ismember(segid,segid(caind))); % sometimes the insertion code is blank
 wgt=mass(resinds); wgt=wgt/sum(wgt);
% length(resinds)
 xpdbcg(caind) = wgt'*xpdb(resinds);
 ypdbcg(caind) = wgt'*ypdb(resinds);
 zpdbcg(caind) = wgt'*zpdb(resinds);
%
 xcg(caind,:) = wgt'*xall(resinds,:);
 ycg(caind,:) = wgt'*yall(resinds,:);
 zcg(caind,:) = wgt'*zall(resinds,:);
end
%
% write out light and heavy chain files separately (as in main)
%
catraf=find(typeCA & traf);
carank=find(typeCA & rnk);
inds={catraf, carank} ;

names={'-cg-traf', '-cg-rank'};

%
for i=1:length(inds)
 ind=cell2mat(inds(i));
 nm=char(names(i));
 disp(['Aligning system ',nm,'...']);
%
 [x, y, z, rmsd]=bestfit(xcg,ycg,zcg,xpdbcg,ypdbcg,zpdbcg,mass,ind);
% compute average structure & realign :
 for j=1:5
  xave=mean(x,2);
  yave=mean(y,2);
  zave=mean(z,2);
  [x, y, z, rmsd]=bestfit(x,y,z,xave,yave,zave,mass,ind);
 end
%
 pdbout([name,nm,'.pdb'], xpdbcg, ypdbcg, zpdbcg, [], [], ind)
% write new dcd :
 writedcd([name,nm,'.dcd'], x(ind,:), y(ind,:), z(ind,:));
end

