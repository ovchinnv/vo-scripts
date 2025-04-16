% sample a psilcs group (with replacement) from from a probability space
if (~exist('nsamp'))
 nsamp = 8 ; % number of samples to draw
end
%
outfile='psilcs-def.str';
%
if (~exist('qnoreplace'))
 qnoreplace=1; % default is to draw with replacement ; however, it might yet be useful to draw without replacement for more aa diversity
end
% array of all residues that can be sampled to create a psilcs group (e.g. peptide)
allpsilcs={ 'ALA' 'ARG' 'ASP' 'GLN' 'LEU' 'THR' 'GLU' 'ILE' 'PHE' 'LYS' 'SER' 'VAL' 'MET' 'ASN' 'PRO' 'TYR' 'HSD' 'GLY' 'TRP' 'CYS' };
% Now, define the probabilities with which the above residues will be chosen ; this is far from an easy task; at minimum, I do not want astronomical possibilities
% rather, would prefer to identify classes of resiudes that would approximately cover the space of 'hit' peptides
% also, I now believe that basic or acidic residues are not useful
pdfpsilcs=[  0     0     0     0     2     1     0     0     2     0     1     1     0     0     0     2     1     0     2     0    ];
%pdfpsilcs=[  0     0     1     0     2     1     0     0     2     1     1     0     0     0     0     2     1     0     2     0    ];
%
%=======================
% check pdf for positivity
if any(pdfpsilcs<0)
 error ( ['PDF is negative at position(s): ', num2str(find(pdfpsilcs<0))] );
end
% crea(te cumulative pdf:
cumdist=cumsum(pdfpsilcs); 
% check for nonnegativity
% normalize
if (cumdist(end) < eps )
 error ( ['The sum of probaililties is: ', num2str(cumdist(end)) ] );
end
cumdist=cumdist/cumdist(end);
%
% draw samples :
%
rands=rand(1,nsamp);
for i=1:nsamp
 inds(i)=find(cumdist>rands(i),1);
 if (qnoreplace) % without replacement (does not make complete mathematical sense, but let's provide this functionality anyway)
  dpdf=cumdist(inds(i));
  if (inds(i)>1) ; dpdf=dpdf-cumdist(inds(i)-1) ; end
  cumdist(inds(i):end)=cumdist(inds(i):end)-dpdf; % note that we could end up with all entries 0 if we have drawn everything
  cumdist=cumdist/cumdist(end) ;% renormalize
 end
end
samples=allpsilcs(inds);
%
fprintf('Sampled ')
fprintf('%s ', samples{:});
fprintf('\n');
% write charmm stream file :
%
fout=fopen(outfile,'w');
fprintf(fout,'%s\n%s\n\n','* psilcs residues','*');
for i=1:nsamp
 fprintf(fout,'set seg%d %s\n',i,samples{i});
end
fprintf(fout,'set nseg %d\nreturn\n',nsamp);
%
fclose(fout);
