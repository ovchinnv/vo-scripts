% read and process complete log from V tesselation calculations
% This is tricky because on a 64 bit machine the record markers are 
% written as 4-byte integers; and they refer to the number of intervening
% 4-byte words
%
% new version supports multiple files (concatenated into one log)
%
% this version processes new logs (code version 1.22.10)

intsize=4 ;     %bytes per int
intfmt='int32'; %integer format (essentially, 32/64 bit)

if ~(exist('read'))
 read=1;
end
if (read) 
%%%%%%%%%%%%%%%%%%%%%%
%numrep=32; % number of replicas (optional)
fnames={'voro105.log','voro106.log','voro107.log','voro108.log',...
        'voro109.log','voro110.log','voro111.log','voro112.log', ...
        'voro113.log','voro114.log','voro115.log',...
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfile=length(fnames);
data=zeros(0,1);
% loop over files
for j=1:nfile
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
tbeg=0;
tbeg=tmin+round ( 0. * (tmax-tmin) );
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
% fix diagonals if zero
du=diag(his,1);  du=min(1,du); du=1-du; %ones in place of zeros; elsewhere -- zeros
dl=diag(his,-1); dl=min(1,dl); dl=1-dl;
his2=diag(du,1)+diag(dl,-1) ;
% add matrix with fixes
his=his+his2;
%
%%%%%%%%%%%%% construct rate matrix
r=zeros(n);
for i=1:n
  r(i,:)= his(:,i)./occupancy;
  r(i,i)= - sum(his(i,:))/occupancy(i);
end
%%%%%%%%%%%%%% set c(1)=1
c=zeros(n,1);
c(1)=1;
f=zeros(n,1);
f(:)=f(:)-r(:,1)*c(1);
c(2:end)=r(1:end-1,2:end)\f(1:end-1);
%%%%%%%%%%%%%%%% solve for FE
beta=0.59582;

f=-beta*log(c);
%
%smooth free energy
alpha=[0:length(f)-1]; alpha=alpha/alpha(end);
ds=3;
f=smooth2(alpha,f,ds); f=f-f(1);
c=exp(-f/beta); % smoothed c
%
close;
figure('position',[200,200,600,250]);
subplot(1,2,1);
plot(alpha,f,'r*-');%hold on;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MFPT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% 2: compute Q matrix
%1 : consider c as a pdf, and, thus, normalize to 1;
c=c/sum(c);

%2 calculate the number of hops between milestones (N), and the time spent at a milestone (R)
whereami =zeros(n,1);
when     =tbeg*ones(n,1);
N=zeros(n,n,n);
R=zeros(n,n);
mstone=zeros(1,n)-1; %mstone(1)=2; mstone(n)=n-1;
for l=1:ncross
 id=data(1,l); i=data(2,l); j=data(3,l); k=data(4,l); t=data(5,l)-1; % note: k is either i or j 
 whereami(id)=k;
 m=mstone(id); %(i,m) is the milestone at which replica id is currently sitting
 if m<0 % initial condition
  m=j ; 
 else
  if j~=m  %crossed to a different milestone
   N(i,m,j)=N(i,m,j)+1;  % number of cell crossings from [i,m] to [i,j] (obviously, through cell i); i - cell[replica]; m - old neighbor; j - new
  end
  R(i,m)=R(i,m)+(t-when(id)); %update residence time at milestone m
 end
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
Req=zeros(n,n); % equilibrium residence time at milestone ij
for i=1:n
 Req(i,:)=c(i)*R(i,:)+(c(:).*R(:,i))' ;% Average time spent at milestone [i,j] ; note that Req is symmetric
% Req(i,:)=(c(:)'+c(i)).*(R(i,:)+R(:,i)') ; % (looks like it matters where replica comes from)
 N(i,:,:)=c(i)*N(i,:,:);
end
%
% change units of R to s (need simulation ts, e.g. 1 fs)
Req=Req*2e-15;
% remove 0s from Req (it's OK because in the quotient below, N will be 0 for these entries)
Req(find(Req==0))=-1;
% compute Q (instantaneous transition rate matrix)
q=zeros(n,n,n);
for k=1:n
 q(:,:,k)=squeeze(N(:,:,k))./Req ;% rate of hopping from i,j to i,k
end
q(find(q<0))=0; % just in case a jump occured in the last entry, in which case N>0 but corresponding R=-1 (as initialized) 
% compute mean time of escape from ms i,j
tau=sum(q,3); % some of these are zero, so exclude them from below:
ind=(find(tau~=0));
tau(ind)=1./tau(ind);
p=zeros(n,n,n);
for k=1:n
 p(:,:,k)=q(:,:,k).*tau(:,:); % probability to hit ms i,k after hitting i,j 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% 3: compute mean first passage time to last milestone
% this requires reindexing q into a MxM 2D matrix
% create map:
% number of milestones visited:

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

%populate 2D rate matrix: loop over the relevant entries in Q;
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
semilogy(tn,'k-*');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

