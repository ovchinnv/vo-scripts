% process glycosylation links
% split the links into two groups :
% (1) link between ASN and sugar and (2) links between sugars
% links between sugars define separate sugar segments
% links between protein and sugars will be implemented as patches between protein and sugar segments

if ~isfield(molecule,'Link') ; return ; end
links=molecule.Link;

% remove "unsupported" links that could interfere with the clustering
% in this case, O-glycosylation at mannose
%
iremove=find( ismember({links(:).AtomName2}, ' O3') | ismember({links(:).AtomName2}, ' O') );
if ~isempty(iremove)
 for i=iremove(end:-1:1) % go in reverse, otherwise could move another removal candidate to a valid position
  links(i)=links(end);
  links=links(1:end-1);
 end
end

aname1={links.remove1};
rname1={links.resName1};
chainid1=char({links.chainID1});
resid1=[links.resSeq1];
insertion1=char({links.iCode1});
%
aname2={links.AtomName2};
rname2={links.resName2};
chainid2=char({links.chainID2});
resid2=[links.resSeq2];
insertion2=char({links.iCode2});

% generate sugar segments
plink = ismember(rname1,'ASN') | ismember(rname2,'ASN') | ismember(rname1,'GLN') | ismember(rname2,'GLN') ; % protein links
slink = not ( plink ) ; % sugar links : assume that there are only two types

% group linked residues by connectivity into segments

nlinks=length(links);
% initial number of clusters ; one cluster per link, numbered sequentially :
clusters=cell(nlinks,1);
for i=1:nlinks
 clusters{i}=i;
end
nc=length(clusters);

qchanged=1;
while qchanged
 qchanged=0;
 for i=1:nc
  linds1=(clusters{i}); % links in ith cluster
  for j=i+1:nc
   linds2=(clusters{j}); % links in jth cluster
% check if links have a common atom
   for k=linds1
    r1=resid1(k); % first link atom in ith cluster
    c1=chainid1(k);
    i1=insertion1(k);
    r2=resid2(k); % 2nd link atom in ith cluster
    c2=chainid2(k);
    i2=insertion2(k);
    for l=linds2
     r3=resid1(l); % first link atom in jth cluster
     c3=chainid1(l);
     i3=insertion1(l);
     r4=resid2(l); % 2nd link atom in jth cluster
     c4=chainid2(l);
     i4=insertion2(l);
     if ( (r1==r3 && c1==c3 && i1==i3)||(r1==r4 && c1==c4 && i1==i4)||(r2==r3 && c2==c3 && i2==i3)||(r2==r4 && c2==c4 && i2==i4) )
% absorb cluster j into cluster i ; set cluster j to empty
      clusters{i}=[ (clusters{i}) (clusters{j}) ];
      clusters{j}=[];
      qchanged=1;
     end %if
    end %for
   end %for
  end %for
 end %for
end %while
%
id=0;
% write pdbs corresponding to the clusters, taking care not to include protein atoms
for i=1:nc
 clinks = clusters{i}; %links in this cluster
 if isempty(clinks)
  continue
 end
 id=id+1; % increment segment count
 segid=sprintf('S%-3d',id);
 inds=zeros(length(pdbh),1);
 fpatch=fopen([strtrim(segid),'.str'],'w');
 fprintf(fpatch, '* polysaccharide links\n');
 fprintf(fpatch, '*\n');
 for cl=clinks;
%  cl
  a1=strtrim(aname1(cl));
  r1=resid1(cl);
  c1=chainid1(cl);
  i1=insertion1(cl);
  a2=strtrim(aname2(cl));
  r2=resid2(cl);
  c2=chainid2(cl);
  i2=insertion2(cl);
% code below uses fields from pdbh
  inds1 = ismember(resid(:),r1) & ismember(chainid(:),c1) ;
  inds2 = ismember(resid(:),r2) & ismember(chainid(:),c2) ;
%
  if ~isempty(insertion)
   inds1=inds1 & ismember(insertion(:),i1) ;
   inds2=inds2 & ismember(insertion(:),i2) ;
  end
%
  inds=inds | inds1 | inds2 ;
% write patch (cases for glycans)
  cmd='';
  r1s=num2str(r1);
  r2s=num2str(r2);
%
% decide which patch to write ; cover only links that involve a 'C1'
  if (ismember(a2,'C1')) % glycan link to C1
   if     (ismember(a1,'O1'))
    cmd=['patch 11aa ',segid,' ',r1s,' ',segid,' ',r2s,' setup warn ! sugar link' ];
   elseif (ismember(a1,'O2'))
    cmd=['patch 12aa ',segid,' ',r1s,' ',segid,' ',r2s,' setup warn ! sugar link' ];
   elseif (ismember(a1,'O3'))
    cmd=['patch 13aa ',segid,' ',r1s,' ',segid,' ',r2s,' setup warn ! sugar link' ];
   elseif (ismember(a1,'O4'))
    cmd=['patch 14aa ',segid,' ',r1s,' ',segid,' ',r2s,' setup warn ! sugar link' ];
   elseif (ismember(a1,'O6'))
    cmd=['patch 16ab ',segid,' ',r1s,' ',segid,' ',r2s,' setup warn ! sugar link' ];
   elseif (ismember(a1,'ND2')) % link from ASN to sugar
% to change segid to match the ASN (this hack works if the peptidoglycan link is always at the top) :
%    segid=sprintf('a%-3d',r1);
    seg=ch2seg(c1);
    if ~isempty(seg)
     cmd=['patch NGLA ',seg,' ',r1s,' ',segid,' ',r2s,' setup warn ! peptidoglycan link' ];
    else
     cmd=['! patch NGLA ',seg,' ',r1s,' ',segid,' ',r2s,' setup warn ! peptidoglycan link for an invalid segment id'];
    end
   elseif (ismember(a1,'NE2')) % link from GLN to sugar
% to change segid to match the GLN (this hack works if the peptidoglycan link is always at the top) :
%    segid=sprintf('a%-3d',r1);
    seg=ch2seg(c1);
    if ~isempty(seg)
     cmd=['patch QGLA ',seg,' ',r1s,' ',segid,' ',r2s,' setup warn ! peptidoglycan link' ];
    else
     cmd=['! patch QGLA ',seg,' ',r1s,' ',segid,' ',r2s,' setup warn ! peptidoglycan link for an invalid segment id'];
    end
   end
  elseif (ismember(a1,'C1')) % reverse residue order
   if     (ismember(a2,'O1'))
    cmd=['patch 11aa ',segid,' ',r2s,' ',segid,' ',r1s,' setup warn ! sugar link' ];
   elseif (ismember(a2,'O2'))
    cmd=['patch 12aa ',segid,' ',r2s,' ',segid,' ',r1s,' setup warn ! sugar link' ];
   elseif (ismember(a2,'O3'))
    cmd=['patch 13aa ',segid,' ',r2s,' ',segid,' ',r1s,' setup warn ! sugar link' ];
   elseif (ismember(a2,'O4'))
    cmd=['patch 14aa ',segid,' ',r2s,' ',segid,' ',r1s,' setup warn ! sugar link' ];
   elseif (ismember(a2,'O6'))
    cmd=['patch 16ab ',segid,' ',r2s,' ',segid,' ',r1s,' setup warn ! sugar link' ];
   elseif (ismember(a2,'ND2')) % link from ASN to sugar
% to change segid to match the ASN (this hack works if the peptidoglycan link is always at the top) :
%    segid=sprintf('a%-3d',r1);
    seg=ch2seg(c2);
    if ~isempty(seg)
     cmd=['patch NGLA ',seg,' ',r2s,' ',segid,' ',r1s,' setup warn ! peptidoglycan link' ];
    else
     cmd=['! patch NGLA ',seg,' ',r2s,' ',segid,' ',r1s,' setup warn ! peptidoglycan link for an invalid segment id' ];
    end
   elseif (ismember(a2,'NE2')) % link from GLN to sugar
% to change segid to match the GLN (this hack works if the peptidoglycan link is always at the top) :
%    segid=sprintf('a%-3d',r1);
    seg=ch2seg(c2);
    if ~isempty(seg)
     cmd=['patch QGLA ',seg,' ',r2s,' ',segid,' ',r1s,' setup warn ! peptidoglycan link' ];
    else
     cmd=['! patch QGLA ',seg,' ',r2s,' ',segid,' ',r1s,' setup warn ! peptidoglycan link for an invalid segment id' ];
    end
   end
  end
%
  if ~isempty(cmd)
   fprintf(fpatch, '%s\n', cmd);
  end
 end
%
 fclose(fpatch);
 inds=find(inds);
 for ii=inds'
  pdbh(ii).segID=segid;
 end
 mol3.Model.Atom=pdbh;
 pdbout(mol3,[strtrim(segid),'.pdb'],xpdb,ypdb,zpdb,[],[],inds);
end
%
