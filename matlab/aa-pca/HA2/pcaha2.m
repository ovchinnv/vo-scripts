%make sure sqlite3.mex can be found
addpath('/home/taly/scripts/matlab/sqlite3/mex-sqlite3')
addpath('../')

db='/home/taly/flurepo/flu_strains.db';
if ~exist('tab')
 tab=95 ; % table in the strains database (90, 95, 97-99)
end

if (~exist('qmsa'))
 qmsa=1
end

if (qmsa)
% select reference strain, based on which the head and transmembrane domain will be deleted
 cmd='select name, sseq_ids, sequence from flu_strains where strain_id = ''AHL90133.1''';
 ref=sqlite3(db,cmd);
 seq0=ref.sequence ;
% locate excision points
% head
 aa='MKVKLLV';
 idel0=strfind(ref.sequence,aa)  ;
 aa='IPSIIQSR';
 idel1=strfind(ref.sequence,aa) + length(aa)-1 ;
%
% TM
 aa='GVYQILA';
 jdel0=strfind(ref.sequence,aa) + length(aa) ;
 jdel1=inf;

%
 if ~exist('maxseq')
  maxseq=inf ; % maximum number of sequences to align (set to inf to align all available ones)
 end

 table=['flu_strains',num2str(tab)];
 cmd= ['select rowid, *, length(sequence) as slen from ',table,' where Nunk=0 and sequence not LIKE ''%J%'''] ; % ignore sequences with 'J' because Matlab does not score J's
 strains=sqlite3(db, cmd) ;

 toalign=struct([]);
 nseq=length(strains);
% sadly, there if no way to vector-copy fields
 for i=1:min(nseq,maxseq)
  toalign(i).Header = strains(i).strain_id ; 
  toalign(i).Name = strains(i).name ; 
% remove head and transmembrane domain
  seq=strains(i).sequence ;
  [score,aln]=nwalign(seq0, seq, 'GLOCAL', 0) ;
%  showalignment(aln);
% cut out head and TM using alignment :
  ipos=0;
  ires=0 ;
  while ((ires < idel0) & (ipos<size(aln,2)))
   ipos=ipos+1;
   if (aln(1,ipos)~='-')
    ires=ires+1;
   end
  end
  idel0a=ipos ;

%
  while ((ires < idel1) & (ipos<size(aln,2)))
   ipos=ipos+1;
   if (aln(1,ipos)~='-')
    ires=ires+1;
   end
  end
  idel1a=ipos ;
%
  while ((ires < jdel0) & (ipos<size(aln,2)))
   ipos=ipos+1;
   if (aln(1,ipos)~='-')
    ires=ires+1;
   end
  end
  jdel0a=ipos ;
%
  while ((ires < jdel1) & (ipos<size(aln,2)))
   ipos=ipos+1;
   if (aln(1,ipos)~='-')
    ires=ires+1;
   end
  end
  jdel1a=ipos ;
%
  aln(3,idel0a:idel1a)='-';
  aln(3,jdel0a:jdel1a)='-';
  toalign(i).Sequence = strrep(aln(3,:),'-',''); % create sequence without deleted parts
 end
%
 tic;
 msa=multialign(toalign) ; % multiply align sequences
 t=toc;
 qmsa=0

 save(['pca',num2str(tab),'ha2.mat'],'strains','msa')

end
%
