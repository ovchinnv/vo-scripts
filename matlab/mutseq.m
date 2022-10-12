%
muts={
 'GFRHQN', 'GFRWQN' ; 
 'VALENQH','VALINQH' ;
 'QINGKLN','QINGILN' ;
};
%
if ~exist('read') ; read=1; end
if (read)
 seqs=fastaread('../round-g1/seqs-g4.fa') ;
 gpen=10; % unclear what to choose
 msa=multialign(seqs,'GAPOPEN',gpen,'EXTENDGAP',0) ; % minimize number of gaps, but allow existing ones to be long
 showalignment(msa) ;
 read=0;
end
% find H3
i0=-1;
for i=1:numel(msa)
  inds=strfind(msa(i).Header,'H3N2');
  if ~isempty(inds)
   i0=i
   break
  end
end

assert(i0>0) ; % make sure reference strain is found

newmsa=msa ;
seq0=msa(i0).Sequence ; % use original seq
for mut=muts'
 primer=mut{1} ;
 inds=strfind(seq0,primer) 
 assert(numel(inds)==1) % make sure the sequence exists and occurs once
 idiffs=find( primer - mut{2} ) % locations of different letters
 for ind=inds(:) % over all mathes, although allowing only one above
% check that there are no dashes in the seq
  for i=1:numel(newmsa)
   disp(['Applying primer ',primer,'==>', mut{2}, ' to sequence "',msa(i).Header,'" at position(s) ', num2str(inds)]);
   seq=newmsa(i).Sequence ;
   subseq=seq(ind:ind+numel(primer)-1) ;
   if strfind(subseq,'-') ;
    warning ['matching subsequence has dashes in sequence ',msa(i).Header,':',subseq ];
   end
   for idiff=idiffs(:) ;
    seq(ind+idiff-1)=mut{2}(idiff) ; % make mutation
    newmsa(i).Sequence=seq; % "save"
   end
  end
 end
end

compmsa=[msa;newmsa];
showalignment(compmsa)


addpath('~/scripts/matlab') ;
fid=fopen('seqs-g5.fa','w')
seqwrite(fid,newmsa,'fasta')

