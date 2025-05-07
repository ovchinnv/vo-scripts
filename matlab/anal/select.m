% atom selections
heavy=ismember(strtrim(segid),'HC');
light=ismember(strtrim(segid),'LC');
ab=ismember(strtrim(segid), {'HA1','HA2'});
typeCA=ismember(aname,'CA');
typeC=ismember(aname,'C');
typeN=ismember(aname,'N');
typeO=ismember(aname,'O');
typeH=ismember(aname,'H');
<<<<<<< HEAD
typeOH2=ismember(aname,'OH2');
=======
typeOH2=ismember(aname,'OH2'); % water oxygen TIP3P
>>>>>>> 777f639 (local changes before merge)
backbone = typeCA | typeC | typeN | typeO | typeH ;
noh=element~='H';
