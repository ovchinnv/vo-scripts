#!/bin/bash

# script to convert a PDB file into multiple PDB files with one subunit & one model per file readable by charmm 
# the header is prepended to each coordinate data set
# water molecules are deleted;  the first subscript  indicates the subunit; the second, the model

# Fields: Atom, Atom No, Space, Atom name, Alt Conf indic, Resname, Space, 
#  Chain Ident, Res Seq No, Spaces, x, y, z, Occup, Temp fact, Spaces, Segment ID
#  FIELDWIDTHS=" 6 5 1 4 1 3 1 1 4 1 3 8 8 8 6 6 6 4" 




awk -v dash="_i"\
    -v mdash=""\
    -v mdl=""\
    -v start=0\
    -v nsub=0\
    -v subID=""\
    -v nres=""  '                                                # this is how you declare global variables

BEGIN \
  {
  segname = ARGV[1];
  outname = ARGV[2];
  if (outname == "") outname=segname;

######## determine output file name ################################

  i=match(outname,"\\.");
  if (i == 0) outname=outname ".pdb"

#  i=match(outname, "\\.*[^\\.]*$",seg);    # split off extension 
#                       outfile=substr(outname,1,i-1) "_i" seg[0]; 
# 		       print outfile;
   i=match(segname,"\\.*[^\\.]*$",seg);    # split off extension
   segid=substr(segname,1,i-1)
		       
###################################################################
  
  ARGV[2] = "";
  
  stop=0;
  numres=0; ! initialize
 
  }

{
    if ($1 == "SEQRES") {                                                           # find number of residues
                           if ($3 ~ /[A-Z]/) { 
			                       subu=$3;    
					       
					       nid=split(subID,IDs," ");
					       last=IDs[nid] 

			                       if (subu != last) 
					        {
					            ++nsub
						    subID=subID " " subu
                                                    nres=nres " " $4    # look in 4th column (multiple subunits present)
						}
					      }	
			   else {
					        if (subID == "")
						 {
						   nsub=1
						   #subID=""
						   nres=nres " " $3                      # look in 3rd column
                                                 }
			        }

# debug:     print  nsub, subID, nres
		       }

 if (nsub>1) dash="_"  # if more than 1 subunit present

 if ($1 == "MODEL") { #)||($1 == "ENDMDL")) {
  mdl=$2         # model number preceded by a connector
  mdash="_"
  start=1;
 }      
    
 if ($1 == "ATOM") {                    # only work on atoms

    start=1;  
			   
    if ($4 == "HIS") sub("HIS", "HSD");                         # replace HIS with HSD
    if (($4 == "ILE") && ($3 == "CD1")) {sub("CD1","CD ");}     # replace CD1 with CD in isoleucine

    
######################### determine no of residues for current segment ##############################
    
    if (nsub>1) {
     currentres=$6         # get residue number for current segment
     currentsub=$5
     nid=split(nres,NRs," ");
     ii=index(subID,currentsub)/2;
     numres=NRs[ii]; # number of residues
    }
    
    else {
          currentres=$5
	  currentsub=""
	  numres=1.*nres    # type cast
    }
    resnam=$4  # name of residue
######################################################################################################
    
#    print numres, currentres 

    if (currentres == numres){ # last carboxy group
                       if ($3 == "OXT") sub("OXT","OT1") 
		        else
		       if ($3 == "O") sub("O  ", "OT2")
    		       }	    
    if (currentsub !="") sub(resnam " " currentsub " ",resnam "   ")   # remove subunit info if present		       
    
                       i=match(outname, "[\\.][^\\.]*$",seg);    # split off extension 
                       outfile=toupper(substr(outname,1,i-1) dash currentsub mdash mdl seg[0]); 
#		       print outfile;

     print substr($0,1,72) segid >outfile		 

    } # if atom		      

    if ($1 == "TER") { # terminate file
		       print "END" >outfile;
                      }
		      
    if (start==0) {  # write header to header file

                       i=match(outname, "[\\.][^\\.]*$",seg);    # split off extension 
                       outfile=toupper(substr(outname,1,i-1) ".hdr"); 
#                      print outfile, nsub,k;
		       print $0 > outfile;
		       
	       	} # write header
}

END \
{
# now concatenate header + coordinate files

      nid=split(subID,IDs," ");
      i=match(outname, "[\\.][^\\.]*$",seg);    # split off extension 
      headfile=toupper(substr(outname,1,i-1) ".hdr");
      
       if (mdl=="") m=""         
       else m=1 
       do {
        for (k=1;k<=nsub;k++)       # loop over all subunits
                {  # write to all files
#		    print IDs[k]
		       outfile=toupper(substr(outname,1,i-1) dash IDs[k] mdash m seg[0]); 
                       system("mv " outfile " " outfile "_TMP");
		       system("cat " headfile " " outfile "_TMP >" outfile)
		       system("rm -f " outfile "_TMP")  
	       	} # write
	   ++m
	   }
	while (m<=mdl)   	
      system("rm -f " headfile);
}
    ' $1 $2 # arguments to awk

