function pp=readspl1d(forfile)
 fid=fopen(forfile,'r');
 lines=textscan(fid,'%s',1, 'Delimiter','\n');
 nz=strread(lines{1}{1},'%d',1);
 z=textscan(fid,'%f',nz) ; z=z{:} ;
 textscan(fid, '%s',1, 'Delimiter','\n'); % comment
%
 lines=textscan(fid,'%s',1, 'Delimiter','\n');
% pdeg=strread(lines{1}{1},'%d',1); % reads as int32, which is not convertible !
 pdeg=strread(lines{1}{1},'%s',1); pdeg=str2num(pdeg{1});
%
 pcoef=textscan(fid,'%f',(nz-1)*4) ; pcoef=reshape(pcoef{:},4,[])';
% textscan(fid, '%s',1, 'Delimiter','\n'); % comment
 pp.breaks=z';
 pp.pieces=nz-1;
 pp.order=pdeg+1;
 pp.form='pp'; % essential
 pp.dim=1;
 pp.coefs=pcoef;
end