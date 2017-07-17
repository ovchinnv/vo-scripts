% read dcd trajectories
if (exist('qdcd')) ; if (qdcd==1) ; return ; end ; end
if (exist('dcdstep')) ; dcdstep=max(1,dcdstep) ; else ; dcdstep=1 ; end
%
iframe=0 ;
for fname = dcdnames
 dcdfile = [char(fname)];
 disp(['==> Processing file ',dcdfile]);
 h=read_dcdheader(dcdfile) ;
% loop over frames
 if (exist('maxframes'))
  nframes=fix(maxframes/dcdstep);
 else
  nframes=fix(h.NSET/dcdstep);
 end
 xdcd=zeros(natom,nframes);
 ydcd=zeros(natom,nframes);
 zdcd=zeros(natom,nframes);
%
 for k=1:nframes
  for kk=1:dcdstep
   [xdcd(:,k),ydcd(:,k),zdcd(:,k)] = read_dcdstep(h) ;
  end
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

