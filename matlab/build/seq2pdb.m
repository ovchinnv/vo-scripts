function seq2pdb(seq,ibeg,pdbname,mins,mchainid,msegid,mx,my,mz)
 ntaa;
 mpdb=struct();
 if (~exist('mx')) ; mx=9999 * ones(1,numel(seq)); end
 if (~exist('my')) ; my=mx; end
 if (~exist('mz')) ; mz=mx; end
 ires=ibeg
 if (~exist('mchainid')) ; mchainid=''; end
 if (~exist('msegid')) ; msegid=''; end
 if (length(msegid)<4) ;
  blank4='    ';
  msegid=[msegid,blank4(1:4-length(msegid))];
 end
 if (~exist('mins')) ; ins=''; else ; ins=mins ; end
 ind=1;
 for a1=seq
  mpdb(ind).AtomSerNo=ind;
  mpdb(ind).AtomNameStruct.chemSymbol='C';
  mpdb(ind).AtomNameStruct.remoteInd='A';
  mpdb(ind).AtomNameStruct.branch='';
  mpdb(ind).AtomName='CA';
  mpdb(ind).chainID=mchainid;
  mpdb(ind).resName=aa3(a1);
  mpdb(ind).resSeq=ires;
%  cins=char(mins);
  mpdb(ind).iCode=ins;
  mpdb(ind).X=mx(ind);
  mpdb(ind).Y=my(ind);
  mpdb(ind).Z=mz(ind);
  mpdb(ind).occupancy=0.0;
  mpdb(ind).tempFactor=0.0;
  mpdb(ind).element='';
  mpdb(ind).charge='';
  mpdb(ind).altLoc='';
  mpdb(ind).segID=msegid;
  ind=ind+1; % atom index
  if (isempty(ins)) ; ires=ires+1; end % residue index
  ins=ins+1;
 end
% print this pdb :
 tmpmol.Model.Atom=fixcharmm(mpdb);
 pdbwrite(pdbname, tmpmol);
end

