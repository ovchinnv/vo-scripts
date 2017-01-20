% postprocess confinement calculations
% produce log files compatible with parseFEP plugin in VMD
%
close all;
clear ;


% obtain force constants (same recipe as FER simulations, and similar to confinement paper)
inds=1:0.5:5 ;
tmass=339.398 ;
ks= ( 0.001 * 1.9.^(inds-1) * 2 * pi ).^2 * tmass

nfiles=length(ks);

dfm=[] ; %averages
dfv=[] ; %variances

dskip=500 ; % initial data to throw away (in addition to first iblock below)
%dskip=0;
step=10;
%

temp=298;
kboltz = 1.987191e-3;
bet=1/temp/kboltz ;

foutnamef = 'confine-f.fepout'; % forward file
foutnameb = 'confine-b.fepout'; % backward file
fidf = fopen(foutnamef,'w');
fidb = fopen(foutnameb,'w');
% write file headers
fprintf(fidf,'%s\n','#            STEP                 Elec                            vdW                    dE           dE_avg         Temp             dG');
fprintf(fidf,'%s\n','#                           l             l+dl             l            l+dl         E(l+dl)-E(l)');
% backward
fprintf(fidb,'%s\n','#            STEP                 Elec                            vdW                    dE           dE_avg         Temp             dG');
fprintf(fidb,'%s\n','#                           l             l+dl             l            l+dl         E(l+dl)-E(l)');


fids=[fidb,fidf];
%

for j=1:2 % backward and forward
 dgtot=0; % overall FE change
 di=j*2-3;
 ibeg = 1*(j-1) + (nfiles)*(2-j) ;
 iend = ibeg + di * (nfiles-2) ;

 fid=fids(j) ;

 for i=ibeg:di:iend
  istr=sprintf('%2.1f',inds(i));
  fname=['fer-',istr,'.dat'];
  d=load(fname);
%
% throw away part of sample (equilibration)
  t=d(dskip+1:end,1) - d(dskip,1); % time
  msdx2=d(dskip+1:end,2); % % mean-squared deviation
  de=msdx2 * ( ks(i+di) - ks(i) ) ;
%
  l1=ks(i)/ks(end);
  l2=ks(i+di)/ks(end);
  l1str=sprintf('%5.3f',l1);
  l2str=sprintf('%5.3f',l2);
  fprintf(['Processing lambda window ',l1str,'-',l2str,'\n']);

% now compute
  m=length(de);
  da=zeros(1,m);
  dx=zeros(1,m);
  dxa=zeros(1,m);
%
  da (1)=de(1) ;
  dxa(1)=exp(-bet*de(1)) ;
% sum up
  for k=2:m
   da(k)  = da (k-1) + de(k) ; %eav
   dxa(k) = dxa(k-1) + exp(-bet*de(k)) ; %gav
  end
% now normalize
  da(:) =da(:) ./[1:m]';
  dxa(:)=dxa(:)./[1:m]';
% take log of exponential average
  dg=-log(dxa(:))/bet;
%
%
%%%%%%%%%%%%%% write this window to file %%
% prepare data
%
  dall=zeros(m,9) ;
  dall(:,1)=t ;
  dall(:,2:5)=0 ;
  dall(:,6)=de ;
  dall(:,7)=da ;
  dall(:,8)=temp;
  dall(:,9)=dg ;

  fprintf(fid,'%s\n',['#NEW FEP WINDOW: LAMBDA SET TO ', num2str(l1),' LAMBDA2 ',num2str(l2)]);
% put in bogus equilibration line
  fprintf(fid,'%s\n',['#0 STEPS OF EQUILIBRATION AT LAMBDA ',num2str(l1),' COMPLETED']);
%
  fprintf(fid,'%s\n','#STARTING COLLECTION OF ENSEMBLE AVERAGE');
  for k=step:step:m
   fprintf(fid, '%s','FepEnergy:');
   fprintf(fid, '%7d%16.4f%15.4f%15.4f%15.4f%15.4f%15.4f%15.4f%15.4f', dall(k,:) ) ;
%  'FepEnergy:'    100     -10406.4437    -10395.6386      1205.9142      1204.0558         8.9467         9.6551       189.7151         9.4993
   fprintf(fid,'\n');
  end
  dgtot=dgtot+dall(end,9); % total FE change
%                   #Free energy change for lambda window [ 0 0.05 ] is 7.54379 ; net change until now is    7.54379
  fprintf(fid,'%s',['#Free energy change for lambda window [ ',num2str(l1),' ',num2str(l2),' ] is ']); 
  fprintf(fid,'%7.5f',dall(end,9)) ;
  fprintf(fid,'%s',' ; net change until now is ');
  fprintf(fid,'%7.5f',dgtot) ;
  fprintf(fid,'\n');

  end % over all windows
 fclose(fid);
end % backward and forward parts

