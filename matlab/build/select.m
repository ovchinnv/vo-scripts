% atom selections
heavy=ismember(segid,'HC');
light=ismember(segid,'LC');
typeCA=ismember(aname,'CA');
typeC=ismember(aname,'C');
typeN=ismember(aname,'N');
typeO=ismember(aname,'O');
typeH=ismember(aname,'H');
backbone = typeCA | typeC | typeN | typeO | typeH ;

