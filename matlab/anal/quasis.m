% quasiharmonic analysis
%
names={'3bnc60t', '3bnc60at', '3bnc60glt'};
names=[ names {'ch103t', 'ch103-i3.2t', 'ch103ucat'}];
names=[ names {'pgt121t', '3h109lt', 'gl121t'} ];


temp=298; %K
% entropy arrays
tsclass=zeros(3,2);
tsquant=zeros(3,2);

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
% quasiharmonic analysis:
% inds=[ {find(heavy & typeCA)} {find(light & typeCA)} ];
 inds=[ {find(heavy & backbone)} {find(light & backbone)} ];
% tag='ca';
 tag='bb';

 for i=1:2
% align
  ind=cell2mat(inds(i));
  [x, y, z, rmsd]=bestfit(xall,yall,zall,xpdb,ypdb,zpdb,mass,ind);
% compute average structure & realign :
  for j=1:5
   xave=mean(x,2);
   yave=mean(y,2);
   zave=mean(z,2);
   [x, y, z, rmsd]=bestfit(x,y,z,xave,yave,zave,mass,ind);
  end
  [tsquant(ii,i), tsclass(ii,i)]=calc_quasi(x(ind,:), y(ind,:), z(ind,:), mass(ind), temp)
 end %i
end %ii

save(['entropy-',tag,'.mat'], 'tsquant', 'tsclass');
