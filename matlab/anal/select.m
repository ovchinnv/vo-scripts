% atom selections
heavy=ismember(strtrim(segid),'HC');
light=ismember(strtrim(segid),'LC');
ab=ismember(strtrim(segid), {'HA1','HA2'});
typeCA=ismember(aname,'CA');
typeC=ismember(aname,'C');
typeN=ismember(aname,'N');
typeO=ismember(aname,'O');
typeH=ismember(aname,'H');
typeOH2=ismember(aname,'OH2');
backbone = typeCA | typeC | typeN | typeO | typeH ;

