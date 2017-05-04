%
% preprocessor script to combine energy data suitable for VMD's parseFEP : a forward FE file and a backward FE file
% leaving some notes behind as I figure out what is actually being output
% also, for now assuming that all runs within a window have the same number of steps, which simplifies averaging
%
%temperature :
temp=298 ;
%
dskip=100 ; % how many entries to skip when computing final average
fepdir='.' ;
finbase='3h109l' ;
foutbase='3h109l' ;
%
%forward files
lambda0=0;
lambda1=1;
nwin=32;
dlambda=1/nwin;

bex=1 ; % exponent to stretch near 1

% ll=`awk "function sgn(n){ return (n<0)? -1:1 } ; function abs(n){ return (n>0)? n:-n }  BEGIN {  printf \"%.*f\", $dec,  (sgn(2*$l-1)*(abs(2*$l-1))^$bexp + 1)/2 }"`

% note : IN THIS PARTICULAR CASE
%can only consider electrostatic switching off because vdw iteractions are identical in vac and in model
lambda0=0.5; nwin=16;

lambdaf=[lambda0 : dlambda : lambda1]; %  forward sims
lambdaf=(sign(2*lambdaf-1).*abs(2*lambdaf-1).^bex + 1)/2;
%return
lambdab=lambdaf(end:-1:1);

lambdas=lambdaf ; % choose direction : lambdaf/lambdab
lambdas=lambdab ;

nrun=zeros(nwin); % define number of output files per window :
irun=zeros(nwin); % define initial run index
% here, can specify instead custom nruns for each window

nrun(:)=1;
irun(:)=0;

% loop over windows and combine data
tempfile='_temp_';
format long ;
kboltz = 1.987191e-3;
bet=1/temp/kboltz ;
istep=1;
deinst=6;
deav=7;
dtemp=8;
dgav=9;
dgtot=0; % overall FE change

foutname = [foutbase,num2str(lambdas(1)),'-',num2str(lambdas(end)),'.fepout'] ;
fid = fopen(foutname,'w');
% write file header
fprintf(fid,'%s\n','#            STEP                 Elec                            vdW                    dE           dE_avg         Temp             dG');
fprintf(fid,'%s\n','#                           l             l+dl             l            l+dl         E(l+dl)-E(l)');

for i=1:nwin
 l1=lambdas(i);
 l2=lambdas(i+1);
 l1str=sprintf('%5.3f',l1); % important to get the correct file name
 l2str=sprintf('%5.3f',l2);
 fprintf(['Processing lambda window ',l1str,'-',l2str,'\n']);

 for j=irun(i):irun(i)+nrun(i)-1 ;% loop over multiple runs per window
%  fname=[fepdir,'/',finbase,'-',num2str(j),'-',l1str,'-',l2str,'.fep'] ; 
  fname=[fepdir,'/',finbase,'-',l1str,'-',l2str,'.fep'] ; 
% preprocess to remove non-numerical fields
%  system(['grep "FepEnergy:" ',fname,' | awk ''{ $1="" ; print $0}'' > ',tempfile]);
%  d=load(tempfile) ;
  data=load(fname);
  dlen=length(data);
  d=zeros(dlen, 9);
  d(:,istep)=data(:,1);
  d(:,deinst)=data(:,3)-data(:,2) ; % inst energy change
  d(:,deav)=cumsum(d(:,deinst))./(1:dlen)' ;
  d(:,dtemp)=temp;
% save minimum energy from first run only :
  if (j==irun(i))
   demin=min(d(:,deinst));
  end
  d(:,dgav) = ( cumsum( exp ( -bet* ( d(:,deinst) - demin ) ) )./(1:dlen)') ; % expoential average but NOT dG -- no log or beta
%
  [m,n]=size(d);
  d2=zeros(m,4) ; % for extra calculations
%
  d2(1,1) = d(1,deav) ;
  d2(1,2) = d(1,dgav) ;
%
% inst samples:
  for k=2:m
   d2(k,1)=(k+1)*d(k,deav) - k*d(k-1,deav) ; % I know that this is correct b/c it mimics/follows the inst. energy
   d2(k,2)=(k+1)*d(k,dgav) - k*d(k-1,dgav) ;
  end % over this sample
%
  if (j==irun(i)) % create first file
   dall=d ;
   d2all=d2 ;
  else
   prevstep=dall(end,1);
   d(:,1)=d(:,1) + prevstep ; % add previous time steps
   dall=[dall; d] ; % append otherwise
   d2all=[d2all; d2] ; % append
  end
 end % over files in this window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% recompute running averages
% first, throw away dskip entries
 if (dskip>0) ;  dall(:,1)=dall(:,1)-dall(dskip,1) ; end ; % reset timesteps to zero
 dall=dall(dskip+1:end,:) ;
 d2all=d2all(dskip+1:end,:) ;
% now recompute
 [m,n]=size(dall);
 d2all(1,3)=d2all(1,1) ;
 d2all(1,4)=d2all(1,2) ;
% sum up
 d2all(:,3)=cumsum(d2all(:,1))./[1:m]';
 d2all(:,4)=cumsum(d2all(:,2))./[1:m]';
% for k=2:m
%  d2all(k,3) = d2all(k-1,3) + d2all(k,1) ; %eav
%  d2all(k,4) = d2all(k-1,4) + d2all(k,2) ; %gav
% end
% now normalize
% d2all(:,3)=d2all(:,3)./[1:m]';
% d2all(:,4)=d2all(:,4)./[1:m]';
% insert eav into main table
 dall(:,deav)=d2all(:,3);
% take log and insert dgav also
 dall(:,dgav) = demin - log(d2all(:,4))/bet;
%
%
%%%%%%%%%%%%%% write this window to file %%
 fprintf(fid,'%s\n',['#NEW FEP WINDOW: LAMBDA SET TO ', num2str(l1),' LAMBDA2 ',num2str(l2)]);
% put in bogus equilibration line
 fprintf(fid,'%s\n',['#0 STEPS OF EQUILIBRATION AT LAMBDA ',num2str(l1),' COMPLETED']);
%
 fprintf(fid,'%s\n','#STARTING COLLECTION OF ENSEMBLE AVERAGE');
 for k=1:m
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

return ;





