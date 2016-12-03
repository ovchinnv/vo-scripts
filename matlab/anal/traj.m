% read dcd trajectories
if (exist('qdcd')) ; if (qdcd==1) ; return ; end ; end
%
iframe=0 ;
for fname = dcdnames
 dcdfile = [char(fname)];
 disp(['==> Processing file ',dcdfile]);
 h=read_dcdheader(dcdfile) ;
% loop over frames
 nframes=h.NSET;
 xdcd=zeros(natom,nframes);
 ydcd=zeros(natom,nframes);
 zdcd=zeros(natom,nframes);
%
 for k=1:nframes
  [xdcd(:,k),ydcd(:,k),zdcd(:,k)]=read_dcdstep(h) ;
  iframe=iframe+1;
  if (~mod(iframe,100))
   disp(['==> Read frame #', num2str(iframe)]);
  end
  if (exist('maxframes'))
   if (iframe>maxframes) ; break ; end ;% break if number of frames exceeded
  end
 end
% append to global coords
 if(exist('xall'))
  xall=[xall xdcd];
  yall=[yall ydcd];
  zall=[zall zdcd];
 else
  xall=xdcd; yall=ydcd; zall=zdcd;
 end
%
end % dcd
%
nall=size(xall,2);
qdcd=1;

