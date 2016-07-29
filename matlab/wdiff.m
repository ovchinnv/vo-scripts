% compute diffusion coefficient from double half-harmonic window simulations
% loosely based on rate computation from Voronoi logs.

% user options :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temperature=300;%temperature used to convert probability to free energy (K)
timestep=2;     %simulation step in _femtoseconds_ needed to compute the rate (in seconds)
fefile='wfe.mat';%ASCII file that contains the pre-computed free energy profile
diffile='wdiff.mat';% ASCII file that will contain the computed diffusion coefficient
timefromzero=1; %set to 1 if the time step was reset to zero in sequential logs
guessnrep=1 ;   %set to 1 if want to guess number of replicas from log files, otherwise set below
numrep=64;      %number of replicas; ignored if guessnrep = 1
guessfilesize=1;%use system calls to determine log file size for fast reading (linux only); otherwise set below
datasize=10000000;% size of log data array (for faster reading) ; ignored if guessfilesize=1
% Window crosing log files below : use lognames = { 'name1' 'name2' ... etc } basename='./window'
lognames={...
%        'window.log',...
 };
% for multiple log names
basename='window';
basext='.log';
ibeg=5;
iend=11;
% modifications below this line should not normally be needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~isempty(basename))
  for i=ibeg:iend
   lognames=[lognames {[basename,num2str(i),basext]}];
  end
end
% might need to change values below if integer sizes change
intsize=4 ;     %bytes per int
intfmt='int32'; %integer format (essentially, 32/64 bit)
kboltz=1.987191d-3;        % boltzmann constant
hline=' ================================================================================================\n';
fprintf(hline);
fprintf(' Will compute diffusion constants from window crossing logs by fitting to FPE for overdamped dynamics.\n');
fprintf(hline);
if ~(exist('read'));read=1;end
if (read)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfile=length(lognames);
data=zeros(0,1);
% loop over files
% keep track of time in case the calculation was restarted from zero each time (relevant for timefromzero=1))
toffset = 0;
ifile=1;
%nfile=1;
for j=ifile:nfile
 fname=char(lognames(j));
 fprintf([' Reading Window crossing log file ',fname,'\n'])
%get file size (linux only)
 if (guessfilesize)
  attrib=['ls -s ',fname];
  [i,attrib]=system(attrib);
% s=strread(attrib, '%s','delimiter',fname);
% s=str2num(char(s(1)))*1024/4; % file size in Kbytes * 1024/4 int per kB (since all entries are integer); OS/hardware - specific
  [s,junk]=strread(attrib, '%d%s');
  datasize=s*1024/4; % file size in Kbytes * 1024/4 int per kB (since all entries are integer); OS/hardware - specific
 end
 fid=fopen(fname,'r');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% `allocate' data: this is essential for fast speed
%data=zeros(filesize,1); % use
 d=zeros(datasize,1);
 n=fread(fid,1,'int32')/intsize; % this record delimiter (record size in Fortran binary)
 d(1:n)=fread(fid,n,intfmt);i=n+1;
 n=fread(fid,1,'int32')/intsize; % this record delimiter
 n=fread(fid,1,'int32')/intsize; % next record delimiter

 while (length(n)>0)
  d(i:i+n-1)=fread(fid,n,intfmt);i=i+n;
  n=fread(fid,1,'int32')/intsize; % this record delimiter
  n=fread(fid,1,'int32')/intsize; % next record delimiter
 end
 % trim data
 d=d(1:i-1);
% ad hoc : fix time : here assume that the data is consecutive, just the timestep numbering is offset
 if (timefromzero)
  times=d(3:3:end) + toffset;      % extract time and correct
  d(3:3:end)=times ;               % replace time
  toffset = max(times) ;           % compute new offset
 end
%
 data=[data;d];
end % loop over files
%
data=reshape(data,3,[]);
ntotcross=length(data)/3; % number of total crossing events logged
% determine number of replicas (optional):
if (guessnrep) ; nrep=max(data(1,:)); else ; nrep=numrep ; end % replica number is in first position
d=data; %want to modify (truncate) data below, but keep the original record (in d)
read=0; %do not reread this data if script rerun
else
 fprintf(' Using previouly read logs. Type "read=1" or "clear" to force re-read.\n');
end % read
%%%%%%%%%%%%%%%% to consider a subset of the trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time=d(3,:);
tmin=min(time); tmax=max(time);
%
%plot(time) ; % check to make sure time increasing correctly
%
%select a window of time
tbeg=tmin+round ( 0 * (tmax-tmin) );
tend=round(tmax);
ind=intersect ( find(time>=tbeg), find(time<=tend) );
data=d(:,ind);
ntotcross=length(ind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(['==> Pruning crossing data\n']);

mfp=zeros(nrep,2);
for j=1:nrep
 inds=find( d(1,:) == j ) ; % crossing events for replica j
 p=d(2,inds) ; % position data
 s=d(3,inds) ; % step # data
 n=length(inds); % size of data
% generate pruned data :
 m=zeros(100,1); % resized
 t=zeros(100,1);
 tskip=0; % offset for pruned trajectries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 mind=0 ; % initialize milestone index
 ind=0; % initialize
%
 while ind < n
  ind=ind+1;
  mnew=p(ind) ; % milestone that was crossed
  if ~mnew ; continue ; end % skip entries into box (0s), and find next OOB entry
% otherwise : 
  if ( ~mind || mnew~=m(mind) )
   mind=mind+1 ; % increment milestone counter
   if (length(m)<mind); m=[m ; zeros(mind,1) ]; t=[t ; zeros(mind,1) ]; end ;
   m(mind)=mnew ; % record milestone
   t(mind)=s(ind)-tskip;
  end
  tout=s(ind) ;    % record time out
%
% find index of crossing back into box;
% this is the time index at which walker is at the first milestone (i.e. throw away all time before this)
% however, it might happen that the RW jumps across the box in one step, if the noise/timestep are high.
% in that case, all we can say that the trip took less than one step ; will use 0 for the passage time,
  while ind < n
   ind=ind+1 ;
   mnew=p(ind);
   if (m(mind)==mnew) % somehow the same milestone was logged
    continue
    error('Two adjacent escapes through the same boundary logged...');
   end
% otherwise
   tin=s(ind) ; % record time in
   tskip=tskip + ( tin - tout ) + 1 ; % increment time spent OOB

%   t(mind)=s(ind)-tskip ; % record time of milestone, subtracting time spend OO
   % check if index is not zero i.e. escaped on the other side
   if mnew % new milestone crossed
    mind=mind+1 ;
    if (length(m)<mind); m=[m ; zeros(mind,1) ]; t=[t ; zeros(mind,1) ]; end ;
    m(mind)=mnew ; % record new milestone
    t(mind)=s(ind)-tskip;
    tout=s(ind) ; % time out same as time in
    % go to next inner while iteration because we are STILL OOB, albeit on the other side
   else
    break ; % go back to outer while loop, because we are now inside the box
   end % if
  end % while
 end % while

% remove zero entries
 if ~t(mind) ; mind=mind-1 ; end % ignore last milestone if we did not actually escape from it BC PT cannot be computed

 t=t(1:mind);
 m=m(1:mind);
% compute mfpts
 dt=diff(t);
%%%%%%%%%%%%%% check if there are no samples
 if (length(dt)<2)
  mfp(j,:)=NaN; % will deal with this later
  warning(['Window ',num2str(j),' has no valid samples. Set both corresponding MFPTs to NaN.']);
  continue
 end
%%%%%%%%%%%%%
 mfp(j, (sign(m(1))+3)/2 ) = mean(dt(1:2:end)) ;  % if m(1) =-2 then dt(1) is the time to reach 2 from -2 ; put this into mfp(1); otherwise put into mfp(2)
 mfp(j, (sign(m(2))+3)/2 ) = mean(dt(2:2:end)) ;
%
 mfpe(j, (sign(m(1))+3)/2 ) = var(dt(1:2:end)) / length(dt(1:2:end));  % if m(1) =-2 then dt(1) is the time to reach 2 from -1 ; put this into mfp(1); otherwise put into mfp(2)
 mfpe(j, (sign(m(2))+3)/2 ) = var(dt(2:2:end)) / length(dt(1:2:end));
%return
end % over all replicas
%
% compute diffusion from mfp and free energy data
%
% use data from open free energy file
load(fefile) ;
g=mean(gam,1) ; % average over all samples

for i=1:nrep
% free energy derivative
% use gamma and mfpt to compute diffusion coefficient
 Df=1/mfp(i,1)/g(i)^2 * (  exp ( g(i) ) - g(i) - 1 ) ;
 Db=1/mfp(i,2)/g(i)^2 * (  exp (-g(i) ) + g(i) - 1 ) ;
 D=0.5*(Db+Df); % combine two measurements
% D=Df ;
% D=Db ;
 if (g(i)<0) % choose the estimate from the longest MFPT
  D=Df;
 else
  D=Db;
 end
%%
%% another way is to compute D from Tab+Tba (roundtrips) (to fourth order in |b-a| and second order in g: 
 Trt=(mfp(i,1)+mfp(i,2)) ; % round trip time
 Trte=mfpe(i,1)+mfpe(i,2) ; % round drip variance
 eg=exp(g(i)) ; Drt = (eg + 1./eg - 2)/Trt/(g(i)^2) ;
                eDrt= (Drt/Trt)^2 * Trte ; % standard variance of D
% Drt = 1./(mfp(i,1)+mfp(i,2)) ; % truncated
 D=Drt ; % it is essentially the same as D
%
 kBT=kboltz*Temp ; % also in fefile
 dlen=dwin/(nrep-1) ; % length scale conversion to normalize diffusion constant ; this implicitly scales the reaction path length to 1.
 tlen=timestep*1e-15 ; % time in seconds

 Diff(i)=D * dlen^2 / tlen ;
 Diffe(i)=eDrt * (dlen^2 / tlen) ^ 2 ;

end % over all windows
% plot diffusion coefficient with error bars
%
% deal with NaNs:
indnan=find(isnan(Diff)) ;
indok=setdiff([1:nrep],indnan) ;
for i=indnan
 Diff(i)=interp1(indok,Diff(indok),i);
 Diffe(i)=interp1(indok,Diffe(indok),i);
end

figure('position',[300 300 600 250]) ;  hold on ; box on;
errorbar(alpha, Diff, sqrt(Diffe),'k.-') ;
%errorbar(alpha, Diff, sqrt(Diffe),'k.-') ;
xlabel('$\alpha$', 'fontsize',14, 'interpreter','latex');
ylabel('$D(\Delta s)^{-2}$', 'interpreter','latex','fontsize',14);
set(gca, 'fontsize',14);
set(gcf,'paperpositionmode','auto');
print(gcf,'-depsc2','wdiff.eps');
% save fe profile and diffusion coeff. profile
%
save(diffile, 'fe','Diff','Diffe');



