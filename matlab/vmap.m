% create replica map from voronoi log
% this version processes new logs (code version 1.22.10)

intsize=4 ;     %bytes per int
intfmt='int32'; %integer format (essentially, 32/64 bit)

if ~(exist('read'))
 read=1;
end
if (read) 
%%%%%%%%%%%%%%%%%%%%%%
%numrep=32; % number of replicas (optional)
fnames={'voro1.log'};%,'voro18.log','voro19.log','voro20.log','voro21.log'};
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
%end % read
%%%%%%%%%%%%%%%% consider a subset of the trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time=d(5,:);
tmin=min(time); tmax=max(time);
%
tbeg=0;
%tbeg=tmin+round ( 0. * (tmax-tmin) );
tend=round(tmax);
ind=intersect ( find(time>=tbeg), find(time<=tend) );
data=d(:,ind);
ncross=length(ind);
%%%%%%%%%%%%% 1: compute free energy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% reconstruct his matrix
n=numrep;
occupancy=zeros(n);
whereami =zeros(n,1);
when     =tbeg*ones(n,1);
rmap     =zeros(n,tend); 
%
for l=1:ncross
 id=data(1,l); 
 i=data(2,l);  j=data(3,l); k=data(4,l); t=data(5,l)-1;
 rmap(id,when(id)+1:t)=i;
 occupancy(id,i)=occupancy(i)+(t-when(id));
 when(id)=t;
 whereami(id)=k;
end
% approximately, add the remaining steps (approximate because we don't know the real tmax)
for id=1:n
 occupancy(id,whereami(id))=occupancy(id,whereami(id))+tmax-when(id);
 rmap(id,when(id)+1:tmax)=whereami(id);
end
%
end % read
%
% map2 should be the same as whereami'
map2=load('voro1.map');
map2=map2(2,:);

reps=1:numrep;
plot(rmap(reps,1:100:end)','.');


