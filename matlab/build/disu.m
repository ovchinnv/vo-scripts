% create disulide patches based on distance criteria
%
if ~exist('seldisu1'); seldisu1=ones(size(anum)); end
if ~exist('seldisu2'); seldisu2=ones(size(anum)); end
if ~exist('disucut') ; disucut=8.5 ; end
if ~exist('qpsfgen') ; qpsfgen=0; end
% need strtrim here because segids were defined with a trailng space (above)
selcys1=seldisu1&ismember(rname,'CYS')&ismember(aname,'SG');
selcys2=seldisu2&ismember(rname,'CYS')&ismember(aname,'SG');
ic1=find(selcys1);
ic2=find(selcys2);
xc1=xpdb(ic1);yc1=ypdb(ic1);zc1=zpdb(ic1); % to control shape (see below)
xc2=xpdb(ic2);yc2=ypdb(ic2);zc2=zpdb(ic2);

scdist2 = bsxfun('minus', xc2(:), xc1(:)' ).^2 + ...
          bsxfun('minus', yc2(:), yc1(:)' ).^2 + ...
          bsxfun('minus', zc2(:), zc1(:)' ).^2;
idisu = ( scdist2 <= disucut^2 ) & (scdist2 > 1e-4 ) ; % no zero distances to within 1e-4 tol
%
[iidisu,jjdisu]=find(idisu);
iidisu=ic2(iidisu); jjdisu=ic1(jjdisu); % convert to indices
indok=unique( sort([iidisu jjdisu],2), 'rows' ); % sort resids within pairs and remove redundancies
if (qpsfgen);
 fp=fopen('add_disulfides.vmd', 'w');
 fprintf(fp,'#!/bin/vmd\n# disulfide bridge patches\n');
else
 fp=fopen('add_disulfides.str', 'w');
 fprintf(fp,'* disulfide bridge patches\n*\n');
end
for i=1:size(indok,1)
 ii=indok(i,1); jj=indok(i,2);
 if (qpsfgen)
  pcmd=['patch DISU ',strtrim(segid{ii}),':', num2str(resid(ii)),' ',strtrim(segid{jj}),':', num2str(resid(jj)), ' '];
 else
  pcmd=['patch disu ',strtrim(segid{ii}),' ', num2str(resid(ii)),' ',strtrim(segid{jj}),' ', num2str(resid(jj)), ' ', 'setup warn'];
 end
 fprintf(fp,'%s\n',pcmd);
end
fclose(fp);
