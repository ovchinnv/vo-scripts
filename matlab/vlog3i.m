% read and process complete log from V tesselation calculations
% This is tricky because on a 64 bit machine the record markers are 
% written as 4-byte integers; and they refer to the number of intervening
% 4-byte words
%
% new version supports multiple files (concatenated into one log)
%
% this version processes new logs (code version 1.22.10)

maxdist=1;% allow collisions between replicas separated by at most d-1 cells
intsize=4 ;     %bytes per int
intfmt='int32'; %integer format (essentially, 32/64 bit)

if ~(exist('read'))
 read=1;
end
if (read) 
%%%%%%%%%%%%%%%%%%%%%%
%numrep=32; % number of replicas (optional)
%
fnames={
'voro.log_38'                                             ...
'voro.log_39' 'voro.log_40'  'voro.log_41' 'voro.log_42'  ...
'voro.log_43' 'voro.log_44'  'voro.log_45' 'voro.log_46'  ...
'voro.log_47' 'voro.log_48'  'voro.log_49' 'voro.log_50'  ...
'voro.log_51' 'voro.log_52'  'voro.log_53' 'voro.log_54'  ...
'voro.log_55' 'voro.log_56'  'voro.log_57'                ...
};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfile=length(fnames);
data=zeros(0,1);
% loop over files
% keep track of time in case the calculation was restarted from zero each time
toffset = 0;
ifile=11;
%
for j=ifile:nfile
 fname=char(fnames(j));
%get file size (linux only)
 attrib=['ls -s ',fname];
 [i,attrib]=system(attrib);
% s=strread(attrib, '%s','delimiter',fname);
% s=str2num(char(s(1)))*1024/4; % file size in Kbytes * 1024/4 int per kB (since all entries are integer); OS/hardware - specific
 [s,junk]=strread(attrib, '%d%s');
 s=s*1024/4; % file size in Kbytes * 1024/4 int per kB (since all entries are integer); OS/hardware - specific
 fid=fopen(fname,'r');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% `allocate' data: this is essential for fast speed
%data=zeros(100000000,1);
 d=zeros(s,1);
 n=fread(fid,1,'int32')/intsize;
 d(1:n)=fread(fid,n,intfmt);i=n+1;
 n=fread(fid,1,'int32')/intsize; % this record
 n=fread(fid,1,'int32')/intsize; % next record

 while (length(n)>0)
  d(i:i+n-1)=fread(fid,n,intfmt);i=i+n;
  n=fread(fid,1,'int32')/intsize; % this record
  n=fread(fid,1,'int32')/intsize; % next record
 end
 % trim data
 d=d(1:i-1);
% aa : fix time
 times=d(5:5:end) + toffset;      % extract time and correct
 d(5:5:end)=times ;               % replace time
 toffset = max(times) ;           % compute new offset
%
 data=[data;d];
end % loop over files
%
data=reshape(data,5,[]);
ncross=length(data)/5;
% determine number of replicas (optional):
numrep=max(data(1,:));
read=0;
d=data; %want to modify (truncate) data below, but keep the original record (in d)
end % read
%%%%%%%%%%%%%%%% consider a subset of the trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time=d(5,:);
tmin=min(time); tmax=max(time);
%
%plot(time) ; % check to make sure time increasing correctly
%return
%
tbeg=tmin+round ( 0 * (tmax-tmin) );
tend=round(tmax);
ind=intersect ( find(time>=tbeg), find(time<=tend) );
data=d(:,ind);
ncross=length(ind);
%%%%%%%%%%%%% 1: compute free energy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  reconstruct his matrix
n=numrep;
his=zeros(n);
occupancy=zeros(n,1);
whereami =zeros(n,1);
when     =tbeg*ones(n,1);
%
for l=1:ncross
 id=data(1,l); 
 i=data(2,l); j=data(3,l); k=data(4,l); t=data(5,l)-1;
 his(i,j)=his(i,j)+1;
% determine occupancy
 whereami(id)=k;
 occupancy(i)=occupancy(i)+(t-when(id));
 when(id)=t;
end
% approximately, add the remaining steps
for id=1:n
 occupancy(whereami(id))=occupancy(whereami(id))+tmax-when(id);
end
%
%=============== write his matrix to file ==================
fname=['his',num2str(ifile),'.dat'];
fid=fopen(fname,'wt');
for i=1:n
 for j=1:n
  fprintf(fid,'%7d ', his(i,j));
 end
 fprintf(fid,'\n');
end
fclose(fid);
%============================================================
% fix diagonals if zero
du=diag(his,1);  du=min(1,du); du=1-du; %ones in place of zeros; elsewhere -- zeros
dl=diag(his,-1); dl=min(1,dl); dl=1-dl;
his2=diag(du,1)+diag(dl,-1) ;
% add matrix with diagonal fix
his=his + 1 * his2;
%%%%%%%%%%%%% construct rate matrix
occupancy=occupancy/mean(occupancy) ; % for better numerical stability
r=zeros(n);
%
for i=1:n
% compute indices
  ibeg = max(1,i-maxdist) ; iend=min(i+maxdist,n) ; inds=[ibeg:iend]; 
  r(i,inds)= his(inds,i)./occupancy(inds);
  r(i,i)= - sum(his(i,inds))/occupancy(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% icluding all crossing events:
%for i=1:n
%  r(i,:)= his(:,i)./occupancy;
%  r(i,i)= - sum(his(i,:))/occupancy(i);
%end
%%%%%%%%%%%%%% set c(1)=1
c=zeros(n,1);
c(1)=1;
f=zeros(n,1);
f(:)=f(:)-r(:,1)*c(1);
c(2:end)=r(1:end-1,2:end)\f(1:end-1);
c=abs(c); % badly scaled matrices can erroneously produce small negative probabilities
%%%%%%%%%%%%%%%% solve for FE
beta=0.59582;
%
f=-beta*log(c);

alpha=[0:length(f)-1]; alpha=alpha/alpha(end);
dfilter=1; % to make a smoother profile set dfilter > 1
fs=smooth2(alpha,f,dfilter);
close;
figure('position',[100,100,900,350]);
subplot(1,2,1);
lw=1.2;
ms=12;
plot(alpha,fs,'r.-', 'linewidth',lw,'markersize',ms);
xlim([0 1]);
%ylim([-15 22]);
box on;
ylabel('\it F(\alpha) (kcal/mol)', 'fontsize',14);
xlabel('\it \alpha', 'fontsize',14);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MFPT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% 2: compute Q matrix
%1 : consider c as a pdf, and, thus, normalize to 1;
c=c/sum(c);

%2 calculate the number of hops between milestones (N), and the time spent at a milestone (R)
maxdist=2; % makes no significant difference
whereami =zeros(n,1);
when     =tbeg*ones(n,1);
N=zeros(n,n,n);
R=zeros(n,n);
mstone=zeros(1,n)-1; %mstone(1)=2; mstone(n)=n-1;
for l=1:ncross
 id=data(1,l); i=data(2,l); j=data(3,l); k=data(4,l); t=data(5,l)-1; % note: k is either i or j 
 whereami(id)=k;
 m=mstone(id); %(i,m) is the milestone at which replica id is currently sitting
%
 if (abs(i-j)<=maxdist)       % crossed to a different milestone that is within the maximum allowed distance 
%                             % note that the above d-restriction will be very problematic when crossings are actually allowed !)
  if m<0 % initial condition
   m=j;
  else
   if j~=m  %crossed to a different milestone
    N(i,m,j)=N(i,m,j)+1;  % number of cell crossings from [i,m] to [i,j] (obviously, through cell i); i - cell[replica]; m - old neighbor; j - new
   end
  end
 else
  j=m ; %if ( i-j > maxdist ) ignore the cross attempt and keep the current milestone
 end
 if (m>=0)
  R(i,m)=R(i,m)+(t-when(id)); %update residence time at milestone m
 end
%
 when(id)=t; % update current time for replica id
 if (k==i);  % stayed in the same cell
  mstone(id)=j; 
 else ;      % crossed (swap i,j) 
  mstone(id)=i;
 end %record the current milestone; 
end
% approximately, add the remaining steps
for id=1:n
 R(whereami(id),mstone(id))=R(whereami(id),mstone(id))+tmax-when(id);
end
%
% `correct' N, R by equilibrium distribution (c)
% in the calculation below, the milestones are indexed : ij as the boundary between cells i and j ; obviously, ji is the same milestone,
% which is why the (equilibrium) residence times matrix Req(i,j) should be symmetric
% Note that this explanation does not apply to N(i,j,k), which counts transitions from milestone ij to ik; for example N(j,i,k) counts
% transitions from milestone ij (eqv. ji) to ik (not the same as jk)
%
Req=zeros(n,n);   % equilibrium residence time at milestone ij
Neq=zeros(n,n,n);
for i=1:n
 for j=1:n
% could do vectorially, but clearer this way
  Req(i,j)=c(i)*R(i,j) + c(j)*R(j,i) ;% Average time spent at milestone [i,j] ; note that Req is symmetric
  Neq(i,j,:)=c(i)*N(i,j,:) ;
 end
end
%
% change units of R to s (need simulation ts, e.g. 1 fs)
Req=Req*2e-15; % 2fs
% remove 0s from Req (it's OK because in the quotient below, N will be 0 for these entries)
zind = find(Req>0);
iflag=min(min( Req(zind) )); % make negative entries of the same magnitude as the positive ones
Req(find(Req==0)) = -iflag;
% compute Q (instantaneous transition rate matrix)
q=zeros(n,n,n);
for k=1:n
 q(:,:,k)=squeeze(Neq(:,:,k))./Req ; % rate of hopping from i,j to i,k (normalized by the residence time)
end
q(find(q<0)) = 0 ; % just in case a jump occured in the last entry, in which case N>0 but corresponding R<0 (as initialized) 
% compute mean time of escape from ms i,j
tau=sum(q,3); % some of these are zero, so exclude them from below:
ind=(find(tau>0));
tau(ind) = 1. / tau(ind);
p=zeros(n,n,n);
for k=1:n
 p(:,:,k)=q(:,:,k).*tau(:,:); % probability to hit ms i,k after hitting i,j 
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% 3: compute mean first passage time to last milestone
% this requires reindexing q into a MxM 2D matrix
% create map:
% number of milestones visited:
%
%[mi mj]=find(Req>0); % milestones visited (based on nonzero residence times)
% find unique
%x=mi>mj; % true if i>j
% exclude milestone that corresponds to the target milestone, i.e. to which the mean first passage time corresponds ; also exclude [0 0]
%mstones2=setdiff([mi.*x mj.*x], [n n-1; 0 0], 'rows'); % double index
%mi=mstones2(:,1); mj=mstones2(:,2) ; % redefine double indices
%milenum=length(mi)  % number of milestones
% reverse lookup
%mstones1=zeros(n,n); for i=1:milenum; mstones1(mi(i), mj(i))=i; end; % single index (so that now :   mstones2(mstones1(mi(i),mj(i)),:)=[mi(i) mj(i)];  )
%
% find relevant entries in q
qq=zeros(size(q));
%for i=1:n; qq(i,:,:)=squeeze(q(i,:,:))+squeeze(q(i,:,:))'; end % symmetrize q so we count edge pairs for which at least one rate is positive  
%for i=1:n; qq(i,:,:)=min(squeeze(q(i,:,:)),squeeze(q(i,:,:))'); end % antisymmetrize q so we count edge pairs for which BOTH rates are positive  
% NOTE that antysymmetrization does not make a difference here for cases in which matrices are nearly singular (or, perhaps, at all)
ind=find(q>0)-1; % nonzero indices of q; subtract 1 to start at 0
% recreate index triplets
kk=floor(ind/n/n);
jj=floor((ind-kk*n*n)/n);
ii=ind-kk*n*n-jj*n;
% indices should start from 1:
ii=ii+1; jj=jj+1; kk=kk+1;
ind3=[ii jj kk]; % these are the triplets for which q is nonzero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find number of relevant milestones 
mstones2=unique([ii jj ; ii kk; jj ii; kk ii],'rows');
mi=mstones2(:,1); mj=mstones2(:,2) ; % define indices
x=mi>mj; % true if i>j
% throw out redundancy; exclude [0 0]
mstones2=setdiff([mi.*x mj.*x], [0 0], 'rows'); % double index
mi=mstones2(:,1); mj=mstones2(:,2) ; % redefine double indices
milenum=length(mi);  % number of milestones
% reverse lookup: single index (so that:   mstones2(mstones1(mi(i),mj(i)),:)=mstones2(i,:) ;  )
mstones1=zeros(n,n); for i=1:milenum; mstones1(mi(i), mj(i))=i; mstones1(mj(i), mi(i))=i; end; % note that mstones1 is symmetric
%
% populate 2D rate matrix: loop over the relevant entries in Q;
%
qq=zeros(milenum,milenum);
for l=1:length(ind)
 i=ii(l); j=jj(l); k=kk(l);
% if ( mstones1(i,j)>0 & mstones1(i,k)>0 ); %this will exclude the target milestone
  qq(mstones1(i,j),mstones1(i,k))=q(i,j,k);
% end
end
% create diagonal ( q_{ii} = -\sum_{j\neq 1} q_{ij} )
qq=qq-diag(sum(qq,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pick target milestone:
% n,n-1
ind=mstones1(n,n-1);
%ind=mstones1(1,2);
% exclude corresponding row + column
qqq=[qq(1:ind-1,1:ind-1) qq(1:ind-1,ind+1:end) ; qq(ind+1:end,1:ind-1) qq(ind+1:end,ind+1:end) ];
% invert matrix
t=-qqq\ones(milenum-1,1);
t=abs(t); % sometimes ill-conditioned matrices produce negative values
% rehash t because we `deleted' a milestone above
t=[t(1:ind-1);0;t(ind:end)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print mean first passage times from main milestones
mstone=zeros(n-2,1);
tn=zeros(n-2,1);
for i=1:n-2
% it is possible that not all of the main milestones are visited; therefore:
 m=mstones1(i, i+1);
 if(m>0); tn(i)=t(m); end
end
subplot(1,2,2);
semilogy(alpha(2:end-1), tn,'k.-', 'linewidth', lw, 'markersize', ms);
box on;
ylabel('\it MFPT(\alpha) (s)', 'fontsize',14);
xlabel('\it \alpha', 'fontsize',14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gcf, 'paperpositionmode', 'auto');
print(gcf, '-dpsc', 'fe-mfpt-2.eps');
print(gcf, '-djpeg100', 'fe-mfpt-2.jpg');
