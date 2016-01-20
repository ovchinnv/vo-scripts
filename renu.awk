# renumber residues

BEGIN {resnum = 0;
       blank="            ";}
{

if ($1 == "ATOM") { 
                    sub($5, substr(blank,1,length($5)-length(++resnum ""))""resnum)
#                    print $5, substr(blank,1,length($5)-length(resnum ""))""resnum 
		   }
print $0
}
