% write coordinate file
%
%
%   fmt='(I10,10X,2(2X,A8),3F20.10,2X,A8,2X,A8,F20.10)'
%  atomid(i), resnum, resname(i), aname(i), r(1:3,i), segid(i), resid(i), wgt(i)
fid=fopen([corname,'.cor'],'w');
fprintf(fid,'%s\n','* CUSTOM CHARMM COORDINATE FILE FOR CG SIMULATION OF CLPX.');
fprintf(fid,'%s\n','* V. OVCHINNIKOV, HARVARD 2012')
d=date;
fprintf(fid,'%s\n',['* DATE : ',d]);
fprintf(fid,'*\n');
fprintf(fid,'%10d',natom); fprintf(fid,'  EXT\n');

for i=1:natom
 fprintf(fid,'%10d',i)                 ;       %atom id
 fprintf(fid,'%10d',i)                 ;       %residue #
 fprintf(fid,'  %-8s',[beadname,num2str(i)]);   %resname
 fprintf(fid,'  %-8s',[beadname,num2str(i)]);     %atomname
 fprintf(fid,'%20.10f',x(i))  ;   %coordinates
 fprintf(fid,'%20.10f',y(i))  ;   %coordinates
 fprintf(fid,'%20.10f',z(i))  ;   %coordinates
 fprintf(fid,'  %-8s',segid)  ;             %segid
 fprintf(fid,'  %-8d',i)      ;            %resid ( 1 bead per residue )
 fprintf(fid,'%20.10f\n',r(i))      ;            %wmain (putting radius there)
end

fclose(fid)
%!cat test.cor
