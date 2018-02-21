% 2/20/18
% generate a nonbonded exclusion list in CHARMM PSF format
% the input consists of the number of atoms in the PSF, 
% and the file name of atoms numbers in each group ; one entire
% group per line.
% The program outputs the nonbonded section for the PSF file,
% whereby each atom in a given group sees none of the atoms
% in the remaining groups ; the goal is to implement RLES 
% (Ovchinnikov & Karplus, in prep, 2018)
% NOTE, that existing nononded sections in the PSF file are not
% read by this program, and so are not taken into account
%

natom=44 ;
nbxfile='groups' ;

if ~exist('natom');
 error "'natom' undefined. Abort."
end

if ~exist('nbxfile');
 error "'nbxfile' undefined. Abort."
end

if ~exist('outfile');
 error "'outfile' undefined. Abort."
end

if ~exist('qext'); % whether to use extended psf format
 qext=0;
end


info = ' ===>';
warning = ' ===> Warning: ';

%read groups
%[glist,d]=textread(nbxfile,'%s%[^\n]', 'bufsize', 100000); % bufsize : max string length
[glist]=textread(nbxfile,'%s', 'delimiter','\n'); % octave format
%
ngrp=length(glist);
% store each group list in a cell
gcell={};
for i=1:ngrp
 gcell=[ gcell ; str2num(char(glist(i)))];
end
% go over all groups and perform some sanity checks
%
for i=1:ngrp
 grp=cell2mat(gcell(i));
 g=sort(grp);
%
 if ~(any(g==grp))
  printf([warning, 'group ', num2str(i), ' atomlist not sorted.\n']);
  gcell(i)=g;
 end
%
 gu=unique(g);
 if (length(gu)<length(g))
  printf([warning, 'group ', num2str(i), ' has repeated indices.\n']);
  gcell(i)=g;
 end
%
 for j=i+1:ngrp
  grp2=cell2mat(gcell(j));
  gint=intersect(grp,grp2);
  if ( length(gint) > 0 )
   error(['groups ',num2str(i),' and ', num2str(j), ' intersect. Abort.']);
  end
 end
%
end
printf('%s %d %s\n',info,ngrp,'groups found...')
%
% now create nonbonded ecxlusion list
%
% (1) list all atoms for which exclusions are defined :
allinds=[];
ginds=[]; % group indices
for i=1:ngrp
 inds=cell2mat(gcell(i));
 allinds=[allinds inds];
 ginds=[ginds, i*ones(1,length(inds))];
end
allinds=sortrows([allinds;ginds]') ; % increasing atom index, followed by group in which it is found
%
if (any(allinds>natom))
 error(['One or more exclusion index exceeds number of atoms. Abort.']);
end

nbxlist=[]; % initialize list
iblo=zeros(1,natom); % indices into nonbond list; iblo(i) points to last nonbonded exclusion of atom 
% loop over all indices
iblo=0;
iatom=1; % current index
for j=1:size(allinds,1)
% advance to index
 ind=allinds(j,1);
 igrp=allinds(j,2);
 for i=iatom:ind-1
  if (i>1) 
   iblo(i)=iblo(i-1);
  else
   iblo(i)=0;
  end
 end
% exclusion list for atom ind :
 exind=setdiff(allinds(:,1),cell2mat(gcell(igrp))) ;
 exind=exind(find(exind>ind)); % list exclusions w/o duplicates in increasing index order
% fprintf([info, ' found ', num2str(length(exind)), ' exclusions for atom ', num2str(ind),'\n']);
% fprintf(['exclusion list for atom',num2str(ind),':', num2str(exind'),'\n'])
 nbxlist=[nbxlist exind'];
 iblo(ind)=length(nbxlist);
 iatom=ind+1;
end
% now finish the list in case the exlusion list does not contain last atom index:
for i=iatom:natom
 iblo(i)=iblo(i-1);
end

% now can write the exclusion list :
fid=fopen(outfile,'w');
if (qext)
  fmt_nbx='%10d%10d%10d%10d%10d%10d%10d%10d\n'; % extended PSF format
else
  fmt_nbx='%8d%8d%8d%8d%8d%8d%8d%8d\n';
end
nnbx=length(nbxlist);
fprintf(fid,'\n%8d',nnbx) ; fprintf(fid,'%s\n',' !NNB')
fprintf(fid,fmt_nbx,nbxlist)
if (nnbx==0 || mod(nnbx,8)>0 ) ; fprintf(fid,'\n') ; end
%print IBLO array
fprintf(fid,fmt_nbx, iblo)
if (natom==0 || mod(natom,8)>0 ) ; fprintf(fid,'\n') ; end
fclose(fid);
