% process fractional coordinate tranformations
%

Sf=zeros(3,3); % Cartesian to fractional
Uf=zeros(3,1);
Sf(:,1)=[molecule.Scale.Sn1];
Sf(:,2)=[molecule.Scale.Sn2];
Sf(:,3)=[molecule.Scale.Sn3];
Uf(:)=[molecule.Scale.Un];
%
Si=inv(Sf); % inverse transformation : fractional to Cartesian
Ui=-Uf;
%
%
fid=fopen(['ucell.str'],'w');
fprintf(fid,'* unit cell transformations \n*\n\n');
fprintf(fid,'set s11 %12.5f\n', Sf(1,1));
fprintf(fid,'set s21 %12.5f\n', Sf(2,1));
fprintf(fid,'set s31 %12.5f\n', Sf(3,1));
fprintf(fid,'set s12 %12.5f\n', Sf(1,2));
fprintf(fid,'set s22 %12.5f\n', Sf(2,2));
fprintf(fid,'set s32 %12.5f\n', Sf(3,2));
fprintf(fid,'set s13 %12.5f\n', Sf(1,3));
fprintf(fid,'set s23 %12.5f\n', Sf(2,3));
fprintf(fid,'set s33 %12.5f\n', Sf(3,3));
%
fprintf(fid,'set u1 %12.5f\n', Uf(1));
fprintf(fid,'set u2 %12.5f\n', Uf(2));
fprintf(fid,'set u3 %12.5f\n', Uf(3));
% inverse transformations
fprintf(fid,'set si11 %12.5f\n', Si(1,1));
fprintf(fid,'set si21 %12.5f\n', Si(2,1));
fprintf(fid,'set si31 %12.5f\n', Si(3,1));
fprintf(fid,'set si12 %12.5f\n', Si(1,2));
fprintf(fid,'set si22 %12.5f\n', Si(2,2));
fprintf(fid,'set si32 %12.5f\n', Si(3,2));
fprintf(fid,'set si13 %12.5f\n', Si(1,3));
fprintf(fid,'set si23 %12.5f\n', Si(2,3));
fprintf(fid,'set si33 %12.5f\n', Si(3,3));
%
fprintf(fid,'set ui1 %12.5f\n', Ui(1));
fprintf(fid,'set ui2 %12.5f\n', Ui(2));
fprintf(fid,'set ui3 %12.5f\n', Ui(3));

fprintf(fid,'return\n');
fclose(fid);
