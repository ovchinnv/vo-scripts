%# Fields: Atom, Atom No, Space, Atom name, Alt Conf indic, Resname, Space, 
%#  Chain Ident, Res Seq No, Spaces, x, y, z, Occup, Temp fact, Spaces, Segment ID
%#  FIELDWIDTHS=" 6 5 1 4 1 3 1 1 4 1 3 8 8 8 6 6 6 4" 
%# matlab also has pdbread:

fname='clpx-3hws-cgla';

cglq=pdbread([fname,'_chr.pdb']); % charge
cglr=pdbread([fname,'_rad.pdb']); % radius
cglm=pdbread([fname,'_mas.pdb']); % mass

atomsq=cglq.Model.Atom;
atomsr=cglr.Model.Atom;
atomsm=cglm.Model.Atom;

natom=length(atomsq);

% extract coordinates from file
x=zeros(1,natom);
y=zeros(1,natom);
z=zeros(1,natom);
q=zeros(1,natom);
r=zeros(1,natom);
m=zeros(1,natom);


for i=1:natom
 atom=atomsq(i);
 atomr=atomsr(i);
 atomm=atomsm(i);

 x(i)=atom.X;
 y(i)=atom.Y;
 z(i)=atom.Z;
 q(i)=atom.tempFactor;
%
 r(i)=atomr.tempFactor;
 m(i)=atomm.tempFactor;

end

tri=DelaunayTri(x',y',z'); % looks like this returns a trirep object -- can run methods as , e.g., tri.edges
edges=tri.edges;
%
% compute edge lenghts
%

dx = x(edges(:,1)) - x(edges(:,2));
dy = y(edges(:,1)) - y(edges(:,2));
dz = z(edges(:,1)) - z(edges(:,2));
ds = sqrt ( dx.^2 + dy.^2 + dz.^2 );
%
% remove edges that are too long -- beyond a certain cutoff
% however, ensure that each node is connected to at least three othre nodes
%
% create connectivity (and bond distance) matrix
bonds=edges;
G=sparse(bonds(:,1),bonds(:,2),1./ds,natom,natom); % note: G is upper triangular since bond list is unique
G=G+G' ; % symmetrize

shortbonds=zeros(0,2);
minbonds =4; % smallest number of allowed bonds between beads
%find the (three) shortest bonds corresponding to every node
clear inds;
for i=1:natom
 [d,inds]=sort(G(i,:),2,'descend'); % note that I am using 1/distance to make the sort work
% the first three indices correspond to the shortest bonds
 for j=1:minbonds
  shortbonds=[shortbonds ; sort([i,inds(j)])];
 end % add bonds
end % over all atoms
shortbonds=sortrows(shortbonds,1); shortbonds=unique(shortbonds,'rows');

% 
%scut=11.; % overfitting ?
scut=8.; % overfitting ?
%scut=5.5; %
inds=find(ds<scut);
bonds=edges(inds,:);
% append cut bonds to the minimal bonds set; sort and remove duplicates
bonds=[bonds ; shortbonds];
bonds=sortrows(bonds,1);
bonds=unique(bonds,'rows');

%bonds=zeros(0,2); % aa to remove all bonds
%
%recompute bond distances (for parameter file)
%
dx = x(bonds(:,1)) - x(bonds(:,2));
dy = y(bonds(:,1)) - y(bonds(:,2));
dz = z(bonds(:,1)) - z(bonds(:,2));
ds = sqrt ( dx.^2 + dy.^2 + dz.^2 );
%

% save data for stiffness fitting (to be done later)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([fname,'.dat'], 'm', 'q', 'r', 'bonds')


% basic output for charmm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write topology file
topname='test';
run top
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write parameter file
topname='test';
% (1) put in dummy default values
kf = 2 * ones ( length(bonds), 1);
% (2) put in optimized values:
fname='param/924/stiff49.dat';
fit=load(fname);
%kf=fit(:,3)*2; % adhoc scale
%kf(:)=mean(kf);
run par
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write psf file
psfname='test';
run psf
%########################## Write compatible coordinate file ##########################
corname=psfname;
run cor

