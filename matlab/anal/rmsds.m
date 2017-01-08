%
names={'3bnc60t', '3bnc60at', '3bnc60glt'};
names=[ names {'ch103t', 'ch103-i3.2t', 'ch103ucat'} ];
names=[ names {'pgt121t', '3h109lt', 'gl121t'} ];

% rmsd arrays
allrmsd=[];
labels={};

for ii=1:length(names)
 name=char(names(ii));

 qpdb=0 ; qpsf=0 ; qdcd=0; clear xall;% to reread

 pdbfile=['../../struc/',name,'_now.pdb']; % in lieu of a "structure" file
% read structure
 struc ;
 select ;
% read trajectories :
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
 traj;
%
 inds=[ {find(heavy & backbone)} {find(light & backbone)} { find(backbone) } ];
 tags=[ {'hc_bb'} {'lc_bb'} {'bb'} ];
%
 for i=1:3
% align
  ind=cell2mat(inds(i));
  tag=char(tags(i));
  [xall, yall, zall, rmsd]=bestfit(xall,yall,zall,xpdb,ypdb,zpdb,mass,ind);

  allrmsd=[allrmsd ; rmsd] ;
  labels=[labels {[name,'-',tag]}];

 end %i
end %ii

save('rmsds.mat', 'allrmsd', 'labels');

