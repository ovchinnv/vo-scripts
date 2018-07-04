% write topology file
%

fid=fopen([topname,'.rtf'],'w');
fprintf(fid,'%s\n','* CUSTOM CHARMM TOPOLOGY FILE FOR CG SIMULATION OF CLPX. V. OVCHINNIKOV, HARVARD 2012')
d=date;
fprintf(fid,'%s\n',['* DATE : ',d]);
fprintf(fid,'*\n');
fprintf(fid,'%5d%5d\n',[99,1]); % version number
fprintf(fid,'\n')
% typical mass entry:
%
%MASS     1 H      1.00800 ! hydrogen which can h-bond to neutral atom
%MASS   215 BRGA3   79.90400  ! TBRE, 1,1,1-dibromoethane

% each bead has a different type
for i=1:natom
 fprintf(fid,'MASS '); 
 j = i + rtf_offset;
 fprintf(fid,'%5d ',j);
 fprintf(fid,'%-7s ',[beadname,num2str(i)]);   %atomtype
 fprintf(fid,'%8.5f',m(i)); % atom mass
 fprintf(fid,'%s\n',char(desc(i)))          ;   % text descriptor
end

% residue topology entries: they are not required for MD, but _are_ needed to append to the topology file
%RESI MG        2.00 ! Magnesium Ion
%GROUP   
%ATOM MG   MG   2.00
%PATCHING FIRST NONE LAST NONE   
for i=1:natom
 fprintf(fid,'\nRESI %-14s',[beadname,num2str(i)]); 
 fprintf(fid,'%14.6f',q(i));                           fprintf(fid,'%s\n',char(desc(i)))          ;   % text descriptor
 fprintf(fid,'GROUP\n');
 fprintf(fid,'ATOM ');
 fprintf(fid,'%-7s ',[beadname,num2str(i)]);
 fprintf(fid,'%-7s ',[beadname,num2str(i)]);
 fprintf(fid,'%14.6f\n',q(i));
 fprintf(fid,'PATCHING FIRST NONE LAST NONE\n');
end



fprintf(fid,'\n')
fprintf(fid,'DEFA FIRS NONE LAST NONE\n\n');
fprintf(fid,'END\n');

fclose(fid)
%!cat test.rtf

