# COMPARES SEQUENCES FROM TWO PDB FILES

BEGIN { out="sequ.dat"
        file1=ARGV[1];
        file2=ARGV[2];
        file3=ARGV[3];
	ARGV[2]="";
	ARGV[3]="";
	resno=0;

### RECOGNIZED AA FOR NUMBERING

        AA["ALA"]=1;AA["ARG"]=1;AA["ASP"]=1;AA["GLN"]=1;AA["LEU"]=1;AA["THR"]=1;AA["GLU"]=1;AA["ILE"]=1;
	AA["PHE"]=1;AA["LYS"]=1;AA["SER"]=1;AA["VAL"]=1;AA["MET"]=1;AA["ASN"]=1;AA["PRO"]=1;AA["TYR"]=1;
	AA["HIS"]=1;AA["HSD"]=1;AA["HSE"]=1;AA["GLY"]=1;AA["TRP"]=1;AA["CYS"]=1;
	
#	print file2
}

{ 
  if ($1 == "SEQRES") 
   {
#### Number the residues #######

start=0
blank="                                                                 ";
numstr=""
for (i=1;i<=NF;i++){
 if ($i in AA) {
      resno +=1;
      if (start == 0) 
       {
        start=index($0,$i)
	numstr=substr(blank,1,start-1+length(file1)+1) numstr
       } # add blank space
      numstr=numstr resno substr("   ",1,4-length(resno ""))  # insert blanks
                }
                   }
 print numstr>out		   

####   NUMBER RESIDUES

   
   thisline=$0 # save line
   
   print file1,$0>out 

   if (file2 != "") {
   
   while (1) {
     oldline=$0
     getline < file2

     if ($0 == oldline) {break;} ! reached end of file2 because two consecutive lines are the same
     if ($1 == "SEQRES") {
                          print file2,$0>out;
                          break
			  }
     			  
             }
   }	     

################## FILE 3

   if (file3 != "") {
   
   while (1) {
     oldline=$0
     getline < file3

     if ($0 == oldline) {break;} ! reached end of file3 because two consecutive lines are the same
     if ($1 == "SEQRES") {
                          print file3,$0>out;
                          break
			  }
     			  
             }
   }	     

  print "*">out

   }

}
