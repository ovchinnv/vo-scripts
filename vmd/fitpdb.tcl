#!/bin/tcl
source tmdAux5.tcl;

# PDB FILE STANDARD:
#   1. |    1 -  6    |   A6    | Record ID (eg ATOM, HETATM)       
#   2. |    7 - 11    |   I5    | Atom serial number                            
#   -  |   12 - 12    |   1X    | Blank                                         
#   3. |   13 - 16    |   A4    | Atom name (eg " CA " , " ND1")   
#   4. |   17 - 17    |   A1    | Alternative location code (if any)            
#   5. |   18 - 20    |   A3    | Standard 3-letter amino acid code for residue 
#   -  |   21 - 21    |   1X    | Blank                                         
#   6. |   22 - 22    |   A1    | Chain identifier code                         
#   7. |   23 - 26    |   I4    | Residue sequence number                       
#   8. |   27 - 27    |   A1    | Insertion code (if any)                       
#   -  |   28 - 30    |   3X    | Blank                                         
#   9. |   31 - 38    |  F8.3   | Atom's x-coordinate                         
#  10. |   39 - 46    |  F8.3   | Atom's y-coordinate                         
#  11. |   47 - 54    |  F8.3   | Atom's z-coordinate                         
#  12. |   55 - 60    |  F6.2   | Occupancy value for atom                      
#  13. |   61 - 66    |  F6.2   | B-value (thermal factor)                    

#set file1 rigor_am/rigor_am_sn_tmd.res
#set file2 rigor_am/rigor_am_sn_tmd1.res
#set file1 rigor_tmd_files/rigor_sn_tmd1.pdb
set file1 conv6pps_hel3_l34+.tmd
set file2 conv6rig_hel3.res
#set file2 rigor_sn_tmd2.pdb

set atoms "";
set targetcoor "";
set weights "";

set fid1 [open $file1 r]
puts "Reading file 1"
foreach line [split [read $fid1] \n] {
    set recordID [string trim [string range $line 0 5]]
    set atomnum  [string trim [string range $line 6 10]]
    set atomname [string trim [string range $line 12 15]]
    set resname  [string trim [string range $line 17 19]]	
    set resnum   [string trim [string range $line 12 15]]
    set xcrd     [string trim [string range $line 30 37]]
    set ycrd     [string trim [string range $line 38 45]]
    set zcrd     [string trim [string range $line 46 53]]
    set occu     [string trim [string range $line 54 59]]
    set beta     [string trim [string range $line 60 65]]
    set segid    [string trim [string range $line 72 75]]
    
    
    if { [string equal $recordID {ATOM}] } { ;#process only ATOM fields
# different selection criteria can be wired below
#      if { $occu == 1.0 }  ;# occupancy
       if { [lsearch -exact "CA N O C" $atomname] > -1 } {;# backbone
        lappend atoms $atomnum
        lappend targetcoor $xcrd $ycrd $zcrd
        lappend weights 1.0
#        puts "$atomnum $segid $resname $atomname $occu"
#        puts "$xcrd $ycrd $zcrd"
#        addatom $atomnum ;#NAMD
      }
    }
}
close $fid1
set natom [llength $atoms];

puts "number of atoms: $natom"

#puts "$targetcoor"
#set targetcoor [translate "1000. 1000. 1000." $targetcoor];
#puts "$targetcoor"
#return;
############################################################
set atoms "";
set targetcoor1 "";

puts "Reading file 2"
set fid1 [open $file2 r]
foreach line [split [read $fid1] \n] {
    set recordID [string trim [string range $line 0 5]]
    set atomnum  [string trim [string range $line 6 10]]
    set atomname [string trim [string range $line 12 15]]
    set resname  [string trim [string range $line 17 19]]	
    set resnum   [string trim [string range $line 12 15]]
    set xcrd     [string trim [string range $line 30 37]]
    set ycrd     [string trim [string range $line 38 45]]
    set zcrd     [string trim [string range $line 46 53]]
    set occu     [string trim [string range $line 54 59]]
    set beta     [string trim [string range $line 60 65]]
    set segid    [string trim [string range $line 72 75]]
    
    
    if { [string equal $recordID {ATOM}] } { ;#process only ATOM fields
# different selection criteria can be wired below
#      if { $occu == 1.0 }  ;# occupancy
       if { [lsearch -exact "CA N O C" $atomname] > -1 } { ;# backbone
        lappend atoms $atomnum
        lappend targetcoor1 $xcrd $ycrd $zcrd
#        puts "$atomnum $segid $resname $atomname $occu"
#        addatom $atomnum ;#NAMD
      }
    }
}
close $fid1
set natom [llength $atoms];
puts "number of atoms: $natom"
#puts "$targetcoor1"
################################ now can align ####

puts "Computing Fit"

set fit [RMSBestFit $targetcoor1 $targetcoor $weights];
puts "[rmsd $targetcoor1 $targetcoor $weights]";
puts "[rmsd $fit $targetcoor $weights]";


