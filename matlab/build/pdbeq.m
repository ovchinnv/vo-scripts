 function yes=pdbeq(atom1, atom2); % check whether atom1 is further in sequence than atom
 % note: only the chain ID, residue number and insertion code are used, atom numbers are ignored
 % this might change in the future, however, atom ordering within a residue is not important
 c1=atom1.chainID;
 r1=atom1.resSeq;
 i1=atom1.iCode;
%
 c2=atom2.chainID;
 r2=atom2.resSeq;
 i2=atom2.iCode;
%
% NOTE : ie either i1 or i2 are empty, ANY operation usign them, including the relational > below
% produces the empty set; this means that "yes" will often be the empty set;
% to avoid this problem, I convert empties to space, and use isspace in place of isempty;
% this works because relational operations using spaces are valid.
 if isempty(i1) ;i1=' '; end
 if isempty(i2) ;i2=' '; end

 no=(c1~=c2) | ...
       (c1==c2) & ( (r1~=r2) | ...
                    (r1==r2) & ( ~isspace(i1) & ( isspace(i2) | ( ~isspace(i2) & i1~=i2 ) ) ) ...
                  );
%c1,c2,c1>c2
%r1,r2,r1>r2
%i1,i2,i1>i2
%isempty(i1)
%isempty(i2)
yes=~no ;
