% process tempering files

file='gdata.dat';

% extract grid
system(['gdata ',restart_file,' > ',file]);

d=load(file);

bet=d(:,1);
esum=d(:,2);
eesum=d(:,3);
eavg=d(:,4);
eeavg=d(:,5);
wgt=d(:,6);
nsampl=d(:,7);

kb=1.987191e-3 ;
temp=1./(kb*bet);

% check if we are to remove earlier stats from another restart file ( must precede the main one in sequence )
if (exist('restart_subtract'))
 system(['gdata ',restart_subtract,' > ',file]);
 d=load(file);
 bet2=d(:,1);
 if ( max(abs(bet2-bet)) > 1e-7 )
  error('Restart file temperature grids do not match');
  return ;
 end

 esum   = esum-d(:,2);
 eesum  = eesum - d(:,3);
 wgt    = wgt - d(:,6);
 nsampl = nsampl - d(:,7);
% recompute averages
 eavg   = esum./wgt ;
 eeavg  = eesum./wgt ;
end
%
%
return ;

% check that averages are correct
fh=99;
figure(fh); hold on;
plot(eavg-esum./wgt,'r'); % ok to double precision
plot(eeavg-eesum./wgt,'g'); % ok to  single precision

