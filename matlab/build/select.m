% atom selections
%heavy=ismember(segid,'HC');
%light=ismember(segid,'LC');
typeCA=ismember(aname,'CA');
typeC=ismember(aname,'C');
typeN=ismember(aname,'N');
typeO=ismember(aname,'O');
typeOH2=ismember(aname,'OH2');
typeH=ismember(aname,'H');
typeOH2=ismember(aname,'OH2'); % TIP3P water oxygen
backbone = typeCA | typeC | typeN | typeO | typeH ;

% atomic element (crude) :
element=cellfun(@(s) s(1), aname);
%
% now populate atomic number array :
%
proton=find(element=='H');
carbon=find(element=='C');
nitrogen=find(element=='N');
oxygen=find(element=='O');
phosphorus=find(element=='P');
sulphur=find(element=='S');
sulfur=sulphur;
heavy=element~='H';
noh=heavy;
