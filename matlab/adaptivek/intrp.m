% interpolate restart files
% currently, semi-manual, to be improved later
restart_file = 'adaptive.restart129.txt';
interp_file = 'adaptive.restart.interp.txt';
% optional file to remove earlier statistics ; comment to disable
%restart_subtract = '.txt';
format long;
gdata;
close all;

oxmin=0.1 ; % old lower limit
oxmax=10  ; % old higher limit
onp=500   ; % old number of points
cg=0.9 ; %damping constant ; must be less than 1
%--------
nxmin=0.3 ; % new lower limit
nxmax=30  ; % new higher limit
nnp=500   ; % new number of points
ninisampl=250; %initial # samples for interpolated points

% recreate old grid and compare to file :
oxmap=linspace(log(oxmin), log(oxmax), onp)';
err = max ( abs ( oxmap - xmap  ) );
fprintf ([ ' Recreated initial grid to precision ', num2str(err)])
if (err<1e-13)
 fprintf('...PASSED\n');
else
 fprintf('...FAILED\n');
end

% create new grid :

nxmap=linspace(log(nxmin), log(nxmax), nnp)';

if (0) % optional check
% recompute weights using cg constant and nsamp :
% brute force, but maybe will look for analytical solution later
 owgts=zeros(1,onp);
 chkinds=1:100:onp; % spot check ...
 for i=chkinds
  fprintf([' Processing gridpoint ', num2str(i),'\n']);
  w=0;
  for j=1:nsampl(i)
   g=1.0-cg/j;
   w=1+g*w;
  end
  owgts(i)=w;
 end
 err = max ( abs ( owgts(chkinds) - wgt(chkinds)  ) )
 fprintf ([ ' Weights error ', num2str(err)])
end

% interpolate data ; use nn to avoid recomputing weights for all samples
nsavg=interp1(xmap, savg, nxmap, 'nearest','extrap') ;
nssavg=interp1(xmap, ssavg, nxmap, 'nearest','extrap') ;
nnsampl=interp1(xmap, nsampl, nxmap, 'nearest','extrap') ;
nwgt=interp1(xmap, wgt, nxmap, 'nearest','extrap') ;

% assign nsamples to new points :
inew=[find(nxmap<log(oxmin)) find(nxmap>log(oxmax)) ];
nnsampl(inew)=ninisampl;
% recompute weights for these indices :
for i=inew(:)'
  fprintf([' Processing gridpoint ', num2str(i),'\n']);
  w=0;
  for j=1:nnsampl(i)
   g=1.0-cg/j;
   w=1+g*w;
  end
  nwgt(i)=w;
end

% compute sums using the new weights :
nssum=nsavg.*nwgt;
nsssum=nssavg.*nwgt;

% plot to compare :
figure(1); clf ;
plot(xmap, savg,'r'); hold on;
plot(nxmap, nsavg,'k'); hold on;

% write gdata portion ; manage headers manually for now.
% # gdata_i (i=1, gridsize) ; for ith bin write : log(kforce (i)), SUM S(i), SUM S(i)^2, AVG S(i), AVG S(i)^2, WGT(i), NUM_SAMPLES(i) 
fout=fopen(interp_file,'w');
for i=1:nnp
 fprintf(fout,[' gdata_',num2str(i),'=(','%20.15f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f',' )\n'] ,nxmap(i), nssum(i), nsssum(i), nsavg(i), nssavg(i), nwgt(i), nnsampl(i));
end

fclose(fout)


