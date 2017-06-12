function pdb=fixcharmm(pdb);

rname    = {pdb.resName}' ;
aname    = {pdb.AtomName}' ;

% change to charmm terminology
resHIS=ismember(rname,'HIS') ;
typeOXT=ismember(aname,'OXT')  ;
typeILECD1=ismember(rname,'ILE') & ismember(aname,'CD1');

for i=find(resHIS)'
 pdb(i).resName='HSD';
end
%
for i=find(typeOXT)'
 pdb(i).AtomName='OT1';
 pdb(i).AtomNameStruct.chemSymbol='O';
 pdb(i).AtomNameStruct.remoteInd='T';
 pdb(i).AtomNameStruct.branch='1';
end
%
for i=find(typeILECD1)'
 pdb(i).AtomName='CD';
 pdb(i).AtomNameStruct.chemSymbol='C';
 pdb(i).AtomNameStruct.remoteInd='D';
 pdb(i).AtomNameStruct.branch='';
end
%remove element
for i=1:length(pdb)
 pdb(i).element='';
end
