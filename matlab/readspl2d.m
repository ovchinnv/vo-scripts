function pp=readspl2d(forfile)
 fid=fopen(forfile,'r');
 lines=textscan(fid,'%s',1, 'Delimiter','\n');
 nx=strread(lines{1}{1},'%d',1);
 x=textscan(fid,'%f',nx) ; x=x{:} ;
 textscan(fid, '%s',1, 'Delimiter','\n'); % comment
%
 lines=textscan(fid,'%s',1, 'Delimiter','\n');
 ny=strread(lines{1}{1},'%d',1);
 y=textscan(fid,'%f',ny) ; y=y{:} ;
 textscan(fid, '%s',1, 'Delimiter','\n'); % comment
%
 lines=textscan(fid,'%s',1, 'Delimiter','\n');
% pxdeg=strread(lines{1}{1},'%d',1); % reads as int32, which is not convertible !
 pxdeg=strread(lines{1}{1},'%s',1); pxdeg=str2num(pxdeg{1});
 lines=textscan(fid,'%s',1, 'Delimiter','\n');
 pydeg=strread(lines{1}{1},'%s',1); pydeg=str2num(pydeg{1});
%
 pcoef=textscan(fid,'%f',(nx-1)*4*(ny-1)*4) ; pcoef=reshape(pcoef{:},4,4,nx-1,ny-1);
% textscan(fid, '%s',1, 'Delimiter','\n'); % comment
 pp.breaks(1)={x'};
 pp.breaks(2)={y'};
 pp.pieces(1)={nx-1};
 pp.pieces(2)={ny-1};
 pp.order(1)={pxdeg+1};
 pp.order(2)={pydeg+1};
 pp.form='pp'; % essential
 pp.dim=2;
 pp.coefs=pcoef;
end
