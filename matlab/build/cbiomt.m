% process crystallographic symmetry fields (biomt)
%
r350=molecule.Remark350 ;
% find the beginning of the sequences
for i=1:size(r350,1)
 a=r350(i,:);
 ind=strfind(a,'BIOMT1');
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
for j=i:size(r350,1)
 a=r350(j,:);
 ind=strfind(a,'BIOMT');
 if isempty(ind)
  break;
 end
end

biomt=r350(i:j,:);

[s, nrot, a1, a2, a3, d]=strread(biomt','%s %d %f %f %f %f'); % rotation matrices
nrot=max(nrot);

% rearrange data as rotations and translations
Arot=reshape( [a1,a2,a3]', 3,3,nrot ); % note : transpose is the inverse, which is what we actually need for CHARMM b/c LEFT hand rule
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

