addpath('~/scripts/matlab/build')
% plot structure scores in from the beta column

names={'H9V11jsd', 'H9V31jsd'};
names={'H5V14jun', 'H5V34jun'}
names={'H7V14bsa', 'H7V34bsa'}
names={'H17V14h32', 'H17V34h32'}

odir='scratch_2001-01-01';

styles={'ro',  'gx', 'bv', 'm*', 'ks', 'c^'};
step=10;
figure(1) ; clf ; hold on ; box on ;

i=0;
bmod=0;
emod=4;
if (1)
for name=names
 i=i+1;
% name
 fbase=[odir,'/',name{:},'/ranked_'];
 betca=[];
 for imod=bmod:emod
  fname=[fbase,num2str(imod),'.pdb']
  mol=readpdb(fname);
  pdb=mol.Model.Atom;
  anum= [pdb.AtomSerNo]';
  aname={pdb.AtomName}' ;
  rname={pdb.resName}' ;
  segid={pdb.segID}' ;
  resid=[pdb.resSeq]' ;
  insertion=char({pdb.iCode});
  occu =[pdb.occupancy]' ;
  bet  =[pdb.tempFactor]' ;
  xpdb =[pdb.X]';
  ypdb =[pdb.Y]';
  zpdb =[pdb.Z]';
  chainid=[pdb.chainID];
  element=[pdb.element];
  natom=length(pdb);
%
  select ;
  inds=find(typeCA(:) & chainid(:)=='A');
  betca=[betca  bet(inds) ];

 end
 betcam=mean(betca,2);
 betcae=std(betca,[],2);
 plot(resid(inds(1:step:end)), betcam(1:step:end),[styles{i},'-'])
 errorbar(resid(inds(1:step:end)), betcam(1:step:end), betcae(1:step:end),styles{i}(1))

 betcamm(i)=mean(betcam) ;
 betcaem(i)=sqrt(mean(betcae.^2)) ;

%return
end

end

figure(1) ;
xlabel('resid');
ylabel('Afold2 Score')
legend(names) ;
print(gcf, '-depsc2', [names{1},'-score.eps']) ;

figure(2) ; clf ;

b=bar(betcamm,'w') ; hold on ;
set(gca, 'xticklabel',names)
ylabel('Afold2 score')

errorbar(1:length(betcamm), betcamm, betcaem)

print(gcf, '-depsc2', [names{1},'-ascore.eps']) ;

