%make sure sqlite3.mex can be found
addpath('/home/taly/scripts/matlab/sqlite3/mex-sqlite3')

%clear;
db='/home/taly/flurepo/flu_strains.db'
if ~exist('tab')
 tab=97 ; % table in the strains database (90, 95, 97-99)
end
%
if ~exist('maxseq')
 maxseq=inf ; % maximum number of sequences to align (set to inf to align all available ones)
end

table=['flu_strains',num2str(tab)];
%strains=sqlite3(db, ['select *, len(sequence) as slen, from flu_strains',num2str(tab)],' where nunk=0') ;
%cmd= ['select rowid, *, length(sequence) as slen from ',table,' where Nunk=0'] ; 
cmd= ['select rowid, *, length(sequence) as slen from ',table,' where Nunk=0 and sequence not LIKE ''%J%'''] ; % ignore sequences with 'J' because Matlab does not score J's
strains=sqlite3(db, cmd) ;

if (~exist('qmsa'))
 qmsa=1
end

if (qmsa)
 toalign=struct([]);
 nseq=length(strains);
% sadly, there if no way to vector-copy fields
 for i=1:min(nseq,maxseq)
  toalign(i).Header = strains(i).strain_id ; 
  toalign(i).Name = strains(i).name ; 
  toalign(i).Sequence = strains(i).sequence ; 
 end
%
 tic;
 fprintf('Performing multiple sequence alignment using %d sequences ...\n',i);
 msa=multialign(toalign) ; % multiply align sequences
 t=toc;
 qmsa=0
end
%
if ~exist('qsqladd') ; qsqladd=0 ; end
if (qsqladd)
% add multiple alignment to database
 cmd=['alter table ',table,' add msa text'];
 try
  sqlite3(db, cmd);
 catch
  warning('SQLite error, aborting (table already exists ?)')
  return
 end
 for i=1:nseq
  cmd=['update ', table, ' set msa="',msa(i).Sequence,'" where rowid = ', num2str(strains(i).rowid)]
  sqlite3(db, cmd);
 end
end
