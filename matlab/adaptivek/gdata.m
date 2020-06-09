% process tempering files

file='gdata.dat';

% extract grid
system([gdatapath,'/gdata ',restart_file,' > ',file]);

d=load(file);

xmap=d(:,1); % this is the parameter (e.g. temperature, kforce) transformed (e.g. using beta rather than T, or w=log(kforce) rather than kforce)
ssum=d(:,2);
sssum=d(:,3);
savg=d(:,4);
ssavg=d(:,5);
wgt=d(:,6);
nsampl=d(:,7);

% check if we are to remove earlier stats from another restart file ( must precede the main one in sequence )
if (exist('restart_subtract'))
 system(['gdata ',restart_subtract,' > ',file]);
 d=load(file);
 xmap2=d(:,1);
 if ( max(abs(xmap2-xmap)) > 1e-7 )
  error('Restart file grids do not match');
  return ;
 end

 ssum   = ssum-d(:,2);
 sssum  = sssum - d(:,3);
 wgt    = wgt - d(:,6);
 nsampl = nsampl - d(:,7);
% recompute averages
 savg   = ssum./wgt ;
 ssavg  = sssum./wgt ;
end
%
%
return ;

% check that averages are correct
fh=99;
figure(fh); hold on;
plot(savg-ssum./wgt,'r'); % ok to double precision
plot(ssavg-sssum./wgt,'g'); % ok to  single precision

