%
% preprocessor script to combine FEP data into two files suitable for VMD's parseFEP : a forward FE file and a backward FE file
% leaving some notes behind as I figure out what is actually being output
% also, for now assuming that all runs within a window have the same number of steps, which simplifies averaging
%
%temperature :
temp=298 ;
%
dskip=0 ; % how many entries to skip when computing final average
fepdir='fepdata' ;
%fepdir='.' ;
foutbase='emr-tpp-fer-lo' ;
%
%forward files
lambda0=0;
lambda1=1;
dlambda=0.025;
bex=1.4 ; % exponent to stretch near 1

lambdaf=[lambda0: dlambda:lambda1]; %  forward sims
lambdaf=1-(1-lambdaf).^bex ;
%return
lambdab=lambdaf(end:-1:1);

lambdas=lambdaf ; % choose direction : lambdaf/lambdab
lambdas=lambdab ;

nw=length(lambdas)-1;

nrun=zeros(nw); % define number of output files per window :
irun=zeros(nw); % define initial run index
% here, can specify instaed custom nruns for each window
nrun(:)=16;
irun(:)=61;

% loop over windows and combine data
tempfile='_temp_';
format long ;
kboltz = 1.987191e-3;
bet=1/temp/kboltz ;
eav=7; % column indices
gav=9;
dgtot=0; % overall FE change

foutname = [foutbase,num2str(lambdas(1)),'-',num2str(lambdas(end)),'.fepout'] ;
fid = fopen(foutname,'w');
% write file header
fprintf(fid,'%s\n','#            STEP                 Elec                            vdW                    dE           dE_avg         Temp             dG');
fprintf(fid,'%s\n','#                           l             l+dl             l            l+dl         E(l+dl)-E(l)');

for i=1:nw
 l1=lambdas(i);
 l2=lambdas(i+1);
 l1str=sprintf('%5.3f',l1); % important to get the correct file name
 l2str=sprintf('%5.3f',l2);
 fprintf(['Processing lambda window ',l1str,'-',l2str,'\n']);

 for j=irun(i):irun(i)+nrun(i)-1 ;% loop over multiple runs per window
  fname=[fepdir,'/',num2str(j),'-',l1str,'-',l2str,'.fep'] ; 
% preprocess to remove non-numerical fields
  system(['grep "FepEnergy:" ',fname,' | awk ''{ $1="" ; print $0}'' > ',tempfile]);
  d=load(tempfile) ;
  [m,n]=size(d);
  d2=zeros(m,4) ; % for extra calculations
%
%            STEP                 Elec                            vdW                    dE           dE_avg         Temp             dG
%                           l             l+dl             l            l+dl         E(l+dl)-E(l)
%FepEnergy:    100     -10406.4437    -10395.6386      1205.9142      1204.0558         8.9467         9.6551       189.7151         9.4993
%
% make sure we understand output :
% it looks like the averages in cols 7,9 are cumulative averages
% we need to recover averages over just the interval ( to be able to combine them later )
% do not forget to take log of the exp. average at the end
  d2(1,1) = d(1,eav) ;
% exponentiate dg to get the exponential average:
  d(:,gav)=exp ( -bet * d(:,gav) ) ;
  d2(1,2) = d(1,gav) ;
%
  for k=2:m
   d2(k,1)=(k+1)*d(k,eav) - k*d(k-1,eav) ; % I know that this is correct b/c it mimics/follows the inst. energy
   d2(k,2)=(k+1)*d(k,gav) - k*d(k-1,gav) ;
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
% recompute averages running averages
% first, throw away dskip entries
 if (dskip>0) ;  dall(:,1)=dall(:,1)-dall(dskip,1) ; end ; % reset timesteps to zero
 dall=dall(dskip+1:end,:) ;
 d2all=d2all(dskip+1:end,:) ;
% now recompute
 [m,n]=size(dall);
 d2all(1,3)=d2all(1,1) ;
 d2all(1,4)=d2all(1,2) ;
% sum up
 for k=2:m
  d2all(k,3) = d2all(k-1,3) + d2all(k,1) ; %eav
  d2all(k,4) = d2all(k-1,4) + d2all(k,2) ; %gav
 end
% now normalize
 d2all(:,3)=d2all(:,3)./[1:m]';
 d2all(:,4)=d2all(:,4)./[1:m]';
% insert eav into main table
 dall(:,eav)=d2all(:,3);
% take log and insert gav also
 dall(:,gav)=-log(d2all(:,4))/bet;
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





