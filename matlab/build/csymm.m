% process crystallographic symmetry fields
%
r290=molecule.Remark290 ;
% find the beginning of the sequences
for i=1:size(r290,1)
 a=r290(i,:);
 ind=strfind(a,'SMTRY');
 if ~isempty(ind)
  break;
 end
end
%
if isempty(ind)
 disp(' Could not find symmetry record ... ');
 return
end
%
% also find last symmetry entry
%
for j=i:size(r290,1)
 a=r290(j,:);
 ind=strfind(a,'SMTRY');
 if isempty(ind)
  break;
 end
end

symm=r290(i:j-1,:);

[s, nrot, a1, a2, a3, d]=strread(symm','%s %d %f %f %f %f'); % rotation matrices
nrot=max(nrot);

% rearrange data as rotations and translations
Arot=reshape( [a1,a2,a3]', 3,3,nrot ); % note : transpose os the inverse, which is what we actually need for CHARMM b/c LEFT hand rule
for irot=1:nrot
 Arot(:,:,irot) = Arot(:,:,irot)' ;
end
%
Tr=reshape(d,3,nrot) ;
%
% write transformation scripts for CHARMM
%
for irot=1:nrot
 fid=fopen(['rotate',num2str(irot),'.str'],'w');
 fprintf(fid,'* roto-translation %d \n*\n',irot);
 fprintf(fid,'coor rota sele @move end matr\n');
 fprintf(fid,'%12.5f %12.5f %12.5f\n', Arot(:,:,irot)');
 fprintf(fid,'coor trans sele @move end xdir %12.5f ydir %12.5f zdir %12.5f\n', Tr(:,irot));
 fprintf(fid,'return\n');
 fclose(fid);
end

