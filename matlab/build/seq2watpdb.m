function tmpmol=seq2watpdb(x,y,z,mchainid,msegid,pdbname,qaddh)
% assume TIP3 charmm compatible format
 assert(numel(x)==numel(y));
 assert(numel(x)==numel(z));
 nres=numel(x);
%
 if (nres==0);
  tmpmol=struct();
  fprintf('==> Nothing to do\n');
 end
%
 if (~exist('msegid')) ; msegid=''; end
 if (~exist('mchainid')) ; mchainid=''; end
 if (~exist('qaddh')) ; qaddh=0 ; end % whether to explicitly add hydrogens
 if (qaddh)
% offsets relative to water oxygen
  dhx=0.9572 * [1 0 0 ; ...
                cos(104.52/180*pi) sin(104.52/180*pi) 0];
 end
% matlab output cannot handle large numbers
% for example : residue #s must be <= 9999 (per PDB standard) ; atom ids must be <= 99999
% so, split into multiple files
 maxres=9999;
 nfiles = ceil ( nres / maxres ) ;
 if (nres>maxres)
  fprintf('==> The number of residues in output is %d, which is larger than %d\n', nres, maxres)
  outnames={}; msegids={};
  for i=1:nfiles
   outnames=[outnames {[pdbname,num2str(i)]}];
   msegids=[msegids {[strtrim(msegid),num2str(i)]}];
  end
  fprintf('==> Will split output into %d pdb files named %s .. %s\n',nfiles,outnames{1},outnames{end});
 else
  outnames={pdbname};
  msegids={strtrim(msegid)};
 end
%
 for i=1:nfiles
  mseg=msegids{i};
  if (length(mseg)<4) ;
   blank4='    ';
   msegids{i}=[mseg,blank4(1:4-length(mseg))];
  end
 end
%
 mpdb=struct();
 ins='';
 ires=0; % residue
 ind=0; % atom
 ifile=1;
 for iatom=1:nres
  ires=ires+1; % increment residue index
  ind=ind+1; % increment atom index
  mpdb(ind).AtomSerNo=ind;
  mpdb(ind).AtomNameStruct.chemSymbol='O';
  mpdb(ind).AtomNameStruct.remoteInd='H';
  mpdb(ind).AtomNameStruct.branch='2';
  mpdb(ind).AtomName='OH2';
  mpdb(ind).chainID=mchainid;
  mpdb(ind).resName='TP3';
  mpdb(ind).resSeq=ires; % assume 1 atom per residue
  mpdb(ind).iCode=ins;
  mpdb(ind).X=x(iatom);
  mpdb(ind).Y=y(iatom);
  mpdb(ind).Z=z(iatom);
  mpdb(ind).occupancy=0.0;
  mpdb(ind).tempFactor=0.0;
  mpdb(ind).element='';
  mpdb(ind).charge='';
  mpdb(ind).altLoc='';
  mpdb(ind).segID=msegids{ifile}(1:4);
% check whether to add hydrogens
  if (qaddh)
   for j=1:2
    js=char(48+j);
    ind=ind+1 ;
    mpdb(ind)=mpdb(ind-1);
    mpdb(ind).AtomSerNo=ind;
    mpdb(ind).AtomNameStruct.chemSymbol='H';
    mpdb(ind).AtomNameStruct.remoteInd=js;
    mpdb(ind).AtomNameStruct.branch='';
    mpdb(ind).AtomName=['H',js];
    mpdb(ind).X=x(iatom)+dhx(j,1);
    mpdb(ind).Y=y(iatom)+dhx(j,2);
    mpdb(ind).Z=z(iatom)+dhx(j,3);
   end % for
  end % if
%
  if (ires==maxres)
% write this pdb :
   tmpmol.Model.Atom=fixcharmm(mpdb(1:ind));
   pdbwrite(outnames{ifile}, tmpmol);
   ires=0; % reset residue index
   ind=0; % reset atom index
   ifile=ifile+1; % next file
  end
 end
% write final file if needed :
  if (mod(ires,maxres)>0) % otherwise already wrote above
% write this pdb :
   tmpmol.Model.Atom=fixcharmm(mpdb(1:ind));
   pdbwrite(outnames{ifile}, tmpmol);
  end
