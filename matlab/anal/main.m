% MATLAB routines for trajectory analysis
%
name='3bnc60at';

name='ch103t';
name='ch103ucat';
name='ch103-i3.2t';

name='gl121t';
name='3h109lt';
name='pgt121t';

pdbfile=['../../struc/',name,'_now.pdb']; % in lieu of a "structure" file

struc ;
% define some selections
select ;

% dcd trajectory names :
basename=['../../dcd/',name];
flag='_now';
ext='dcd';
%
inds=1:5 ;
dcdnames={};
for i = inds
 dcdnames=[dcdnames { [basename,'-',num2str(i),flag,'.',ext] } ];
end
%
% read dcds :
traj;
% to compute rmsd :
%rmsd=calc_rmsd(xall,yall,zall,xpdb,ypdb,zpdb,mass,find(typeCA));
% generate aligned dcd structures for further analysis using CA atoms and LC or HC

cahc=find(typeCA & heavy);
calc=find(typeCA & light);
inds={cahc, calc} ;

names={'-ca-hc', '-ca-lc'};

for i=1:length(inds)
 ind=cell2mat(inds(i));
 nm=char(names(i));
%
 [x, y, z, rmsd]=bestfit(xall,yall,zall,xpdb,ypdb,zpdb,mass,ind);
% compute average structure & realign :
 for i=1:5
  xave=mean(x,2);
  yave=mean(y,2);
  zave=mean(z,2);
  [x, y, z, rmsd]=bestfit(x,y,z,xave,yave,zave,mass,ind);
 end
%
 pdbout([name,nm,'.pdb'], xpdb, ypdb, zpdb, [], [], ind)
% write new dcd :
 writedcd([name,nm,'.dcd'], x(ind,:), y(ind,:), z(ind,:));
end
