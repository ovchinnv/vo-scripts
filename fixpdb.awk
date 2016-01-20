# This file is fixpdb.awk.
# Useage awk -f fixpdb.awk [segid=wxyz] [chainID=X]   <pdbfile.in >file.out
#                                       [resname=abc] 
# Extracts segments from pdb files and converts to a format acceptable by charmm.
# In command line can specify up to a four character segid with wxyz, e.g. prot. This 
#  field is ignored by current CHARMM versions, but needed for older versions. 
# Can specify a one character chainID. If is specified on command line, extracts
#  only lines whose character in column 22 matches chainID X. Use to extract specific 
#  subunit from pdb file.
# Instead, can specify a three character resname to select HOH or ligands like ARA.
# If resname is specified, extracts only lines whose resname in columns 18-20 
#  matches resname abc value.
# Writes header line as a remark.
# Ignores all other lines not beginning with ATOM or HETATM.
# If a single coordinate value for an atom is present, takes that. 
# If multiple coordinates are present, signified by A, B, .. in column 17, takes only A.
# If protein and HOH lines are present and protein lacks a chainID, takes the 
#  protein lines only.
# Converts HOH to TIP and adds a 3, making TIP3, HIS to HSD, CD1 to CD_ for ILE, 
#  adds the segid in columns 73-76. Converts OXT or OCT1 to OT1 and OCT2 to OT2.
# Renumbers atoms starting from 1.
# Fields: Atom, Atom No, Space, Atom name, Alt Conf indic, Resname, Space, 
#  Chain Ident, Res Seq No, Spaces, x, y, z, Occup, Temp fact, Spaces, Segment ID

BEGIN {FIELDWIDTHS=" 6 5 1 4 1 3 1 1 4 1 3 8 8 8 6 6 6 4"} 
{
	if ($1 == "HEADER")
		print "REMARK" substr($0, 7, 69)
	if ($1 != "ATOM  " && $1 != "HETATM")
		endif	
	else if ($5 != " " && $5 != "A")
		endif
	else if ($6 == resname || $8 == chainID || ($8 == " " && $1 != "HETATM")) 
	{
		atomno++
		if ($6 == "HOH")
		{	$4 = " OH2"
			$6 = "TIP"
			$7 = "3"
		}
		if ($1 == "HETATM")
			$1 = "ATOM  "
		if ($6 == "HIS")
			$6 = "HSD"
		if ($6 == "ILE" && $4 == " CD1")
			$4 = " CD "
		if ($4 == " OXT" || $4 == "OCT1") 
			$4 = " OT2"
		if ($4 == "OCT2")
			$4 = " OT1"
		printf "%6s",$1
		printf "%5d", atomno
		printf "%1s", " "
		printf "%4s", $4
		printf "%1s", " "
		printf "%3s", $6
		printf "%1s", $7
		printf "%1s", " "
		printf "%4s", $9
		printf "%4s", "    "
		printf "%8s", $12
		printf "%8s", $13
		printf "%8s", $14
		printf "%6s", $15
		printf "%6s", $16
		printf "%6s", "      "
		printf "%4s\n", segid
	}

}
END {printf "%3s\n", "END"}
