% create psf file
%

fid=fopen([psfname,'.psf'],'w');
fprintf(fid,'%s\n','PSF EXT')
fprintf(fid,'\n')
fprintf(fid,'%10d',3) ; fprintf(fid,'%s\n',' !NTITLE')
fprintf(fid,'%s\n','* CUSTOM CHARMM STRUCTURE FILE FOR CG SIMULATION OF CLPX.')
fprintf(fid,'%s\n','* V. OVCHINNIKOV, HARVARD 2012')
d=date;
fprintf(fid,'%s\n',['* DATE : ',d]);
fprintf(fid,'\n')
fprintf(fid,'%10d',natom) ; fprintf(fid,'%s\n',' !NATOM')
% now print atom info
%segid=atomsq(1).segID;
for i=1:natom
 fprintf(fid,'%10d ',i)
 fprintf(fid,'%-8s ',segid)  ;   %segid
 fprintf(fid,'%-8d ',i)                 ;   %resid ( 1 bead per residue )
 fprintf(fid,'%-8s ',[beadname,num2str(i)]);   %resname
 fprintf(fid,'%-8s ',[beadname,num2str(i)]);   %atomname
 j=i+rtf_offset;
 fprintf(fid,'%4d ',j)                  ;   %atom type ( each bead is of a different type)
 fprintf(fid,'%14.6f%14.6f',q(i),m(i))  ;   %charge and mass
 fprintf(fid,'%8d',0)                 ;   % imove array (whether atom is fixed, not used here)
 fprintf(fid,'%s\n',char(desc(i)))          ;   % text descriptor
end
% extended :
%i10,1x,a8,1x ,a8 ,1x,a8,1x,a8,1x   ,a6,1x,2g14.6        ,i8 ! xplor
%i10,1x,a8,1x ,a8 ,1x,a8,1x,a8,1x   ,i4,1x,2g14.6        ,i8 ! charmm
%ii    ,lsegid,lresid,lres ,atype(i),atc_x,cg(i),amass(i),imove(i)
%
fmt_bond='%10d%10d%10d%10d%10d%10d%10d%10d\n';
fmt_angl='%10d%10d%10d%10d%10d%10d%10d%10d%10d\n';
fmt_dihe=fmt_bond;
fmt_impr=fmt_bond;
fmt_nnbx=fmt_bond;
fmt_grp=fmt_angl;
% BONDS
%
nbonds=length(bonds);
fprintf(fid,'\n')
fprintf(fid,'%10d',nbonds) ; fprintf(fid,'%s\n',' !NBOND: bonds')
fprintf(fid,fmt_bond,bonds')
if (nbonds==0 || mod(nbonds,4)>0 ) 
 fprintf(fid,'\n')
end
% ANGLES
nangle=0
fprintf(fid,'\n%10d',nangle) ; fprintf(fid,'%s\n',' !NTHETA: angles')
if (nangle==0 || mod(nangle,3)>0 ) 
 fprintf(fid,'\n')
end
% DIHEDRALS
ndihe=0
fprintf(fid,'\n%10d',ndihe) ; fprintf(fid,'%s\n',' !NPHI: dihedrals')
if (ndihe==0 || mod(ndihe,2)>0 ) 
 fprintf(fid,'\n')
end
% IMPROPERS
nimp=0
fprintf(fid,'\n%10d',nimp) ; fprintf(fid,'%s\n',' !NIMPHI: impropers')
if (ndihe==0 || mod(nimp,2)>0 ) 
 fprintf(fid,'\n')
end
% HBDONORS
ndon=0
fprintf(fid,'\n%10d',ndon) ; fprintf(fid,'%s\n',' !NDON: donors')
if (ndon==0 || mod(ndon,4)>0 ) 
 fprintf(fid,'\n')
end
% HBACCEPT
nacc=0
fprintf(fid,'\n%10d',nacc) ; fprintf(fid,'%s\n',' !NACC: acceptors')
if (nacc==0 || mod(nacc,4)>0 ) 
 fprintf(fid,'\n')
end

% NB EXCLUSIONS
% all pairs in a monomer excluded -- C(n,2) interactions
%
nnbx = natom * (natom-1) / 2 ;
xlist=[];
iblo=zeros(1,natom);
offset=0;
for i=1:natom
 xlist=[xlist, i+1:natom];
 offset=offset+(natom-i);
 iblo(i)=offset;
end
%
fprintf(fid,'\n%10d',nnbx) ; fprintf(fid,'%s\n',' !NNB')
fprintf(fid,fmt_nnbx,xlist)
if (nnbx==0 || mod(nnbx,4)>0 ) ; fprintf(fid,'\n') ; end

%print IBLO array (indexed by atoms) : points to last nonbonded exclusion of atom i
fprintf(fid,fmt_nnbx, iblo)
if (natom==0 || mod(natom,8)>0 ) ; fprintf(fid,'\n') ; end
%
% GROUPS (define groups are composed of one particle each)
% NOTE that above, the pointer to first atom in each group (iblo) starts from 0, not 1; this was a bug
% natom=0 % aa
ngrp=natom; nst2=0;
fprintf(fid, '\n%10d%10d',[ngrp, nst2]) ; fprintf(fid,'%s\n',' !NGRP NST2')
fprintf(fid, fmt_grp, [ 0:natom-1 ; 2*abs(q(1:natom)) ; zeros(1,natom) ] )
if (ngrp==0 || mod(ngrp,3)>0 ) ; fprintf(fid,'\n') ; end
%
% LONE PAIRS
numlp=0; numlph=0
fprintf(fid, '\n%10d%10d',[numlp, numlph]) ; fprintf(fid,'%s\n',' !NUMLP NUMLPH')
%
% CROSS TERMS (CMAP) -- not needed if CMAP not specified
%
%numcross=0;
%fprintf(fid, '\n%10d',numcross) ; fprintf(fid,'%s\n',' !NCRTERM')
%fprintf(fid,'\n')

fclose(fid)
%!cat test.psf
%return
