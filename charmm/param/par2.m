% write parameter file
% individual radii for vdw params
%

fid=fopen([topname,'.prm'],'w');
fprintf(fid,'%s\n','* CUSTOM CHARMM PARAMETER FILE FOR CG SIMULATION OF CLPX. V. OVCHINNIKOV, HARVARD 2012')
d=date;
fprintf(fid,'%s\n',['* DATE : ',d]);
fprintf(fid,'*\n\n');
% typical bond entry:
%BOND
%C    C      450.0       1.38!  FROM B. R. GELIN THESIS AMIDE AND DIPEPTIDES

%CG1N1  NG1T1  1053.00     1.1800 ! ACN, acetonitrile; 3CYP, 3-Cyanopyridine (PYRIDINE pyr-CN) (MP2 by kevo)
fprintf(fid,'BOND\n');
for i=1:length(bonds)
 fprintf(fid,'%-6s ',[beadname,num2str(bonds(i,1))]);
 fprintf(fid,'%-6s ',[beadname,num2str(bonds(i,2))]);
 fprintf(fid,'%14.6f ',kf(i));
 fprintf(fid,'%14.6f\n',ds(i));
end

% FROM: nbonds.doc
%    B) Warning Distance Specifications
%warning-dist::=
%        [ WMIN   real ] [ WRNMXD real ]
%
% WRNMXD - keyword defines a warning cutoff for maximum atom displacement from
%          the last close contact list update (used in EXTEnded)
%

%
%NONBONDED  NBXMOD 5  ATOM CDIEL SHIFT VATOM VDISTANCE VSHIFT -
%     CUTNB 8.0  CTOFNB 7.5  CTONNB 6.5  EPS 1.0  E14FAC 0.4  WMIN 1.5
%
%V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
%
%epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)
%Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j
%
%atom  ignored    epsilon      Rmin/2   ignored   eps,1-4       Rmin/2,1-4
%
%carbons
%C      0.000000  -0.110000     2.000000 ! ALLOW   PEP POL ARO

%fprintf(fid,'\nNONBONDED  NBXMOD 5 ATOM RDIEL SHIFT VATOM VDISTANCE VSHIFT CUTNB 11.0  CTOFNB 10  CTONNB 9  EPS 1.0  E14FAC 0.0  WMIN 1.5');
fprintf(fid,'\nNONBONDED  NBXMOD 5 ATOM RDIEL SHIFT VATOM VDISTANCE VSHIFT CUTNB 12.0  CTOFNB 11  CTONNB 9  EPS 1.0  E14FAC 0.0  WMIN 1.5');

%specific parameters
for i=1:natom
 fprintf(fid,['\n',beadname,num2str(i)]); 
 fprintf(fid,'%12.8f',[0.0 -0.2 1.4*r(i)]) ; %12/2012 ad hoc scaling
end
% generic parameters -- disable is specific used
%fprintf(fid,['\n',beadname,'*']); fprintf(fid,'%12.8f',[0.0 -0.25 mean(r)]) ; does not perform well
%fprintf(fid,['\n',beadname,'*']); fprintf(fid,'%12.8f',[0.0 -0.5 3.1]) ; %12/2012 
%fprintf(fid,['\n',beadname,'*']); fprintf(fid,'%12.8f',[0.0 -0.3 3.1]) ; %12/2012 : smaller interaction, same results as above



fprintf(fid,'\nEND\n');


fclose(fid)

%!cat test.prm
%return
