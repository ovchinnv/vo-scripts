#!/bin/tcl
# Script for performing TMD-TCL in NAMD using TCLForces
# Files: tmd.tcl (this file -- main routine) & tmdAux.tcl (subroutines)
# Author: Victor Ovchinnikov, Harvard University/MIT, June 2007

# 6.28.07; Orientation & forcing are now allowed on different sets of atoms, specified in occupancy & beta, respectively
# 8/14/07: periodically output forces on forced atoms

############################################################################
source tmdAux6.tcl;

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

set TMDatoms {};
set TMDtargetCoor {};
set targetcoor {}; # we need this to be static, therefore global

set TMDforcedFlags {}; # contains atoms for orientation
set TMDorientFlags {}; # contains atoms for forcing
set TMDforcedAtoms 0 ;
set TMDorientAtoms 0;
set TMDforcesOut 0; # file ID of force output file

set fid [open $TMDfile r]
foreach line [split [read $fid] \n] {
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
      if { ( $occu == 1.0 ) || ( $beta == 1.0 ) } { ;# orientation based on occupancy; forcing based on beta
        lappend TMDatoms $atomnum
        lappend TMDtargetCoor $xcrd $ycrd $zcrd
        addatom $atomnum ;
        if { $occu == 1 } {
         lappend TMDorientFlags 1.0;
         incr TMDorientAtoms
        } else {
         lappend TMDorientFlags 0.0;
        } 
        if { $beta == 1 } {
#different criteria:   if { ([lsearch -exact "ATP MG" $segid] > -1) && ([string index $atomname 1] != "H") }
         lappend TMDforcedFlags 1.0 ;# want to orient the whole structure, but apply force only to the ATP !
         incr TMDforcedAtoms;
        } else {
         lappend TMDforcedFlags 0.0;
        } 
       }
    }
}
close $fid;# done with file

set TMDnatom [llength $TMDatoms];
if { ($TMDnatom == 0) || ($TMDforcedFlags == 0) || ($TMDorientFlags == 0) } {
    set TMDon 0
    print "TMD-TCL WARNING: no TMD atoms found, no forcing, or no orientation atoms found!";
} else {
    set TMDon 1
    print "TMD-TCL: $TMDnatom TMD atoms found"
    print "TMD-TCL: orienting structure based on $TMDorientAtoms atoms."; 
    print "TMD-TCL: forcing $TMDforcedAtoms atoms with force constant $TMDk (applied to each atom)"; 
    if { [file exists $TMDlogFile] } { file copy -force -- "$TMDlogFile" "${TMDlogFile}.BAK"; }
    set TMDout [open $TMDlogFile w+];
    puts $TMDout "TMD-TCL: $TMDnatom TMD atoms found"
    puts $TMDout "TMD-TCL: orienting structure based on $TMDorientAtoms atoms."; 
    puts $TMDout "TMD-TCL: forcing $TMDforcedAtoms atoms with force constant $TMDk (applied to each atom)"; 
# do we need to get TargetRMSD from an existing log file?
   if { $TMDtargetRMS == -2 } {
    set TMDoutOld [open $TMDlogFileOld r];
    while { [gets $TMDoutOld line] > -1 } {
      set lastline [split $line "\t"];
     }
     close $TMDoutOld
     set TMDtargetRMS [lindex $lastline 4];
     puts $TMDout "TMD-TCL: Obtained initial target RMS value from file $TMDlogFileOld: $TMDtargetRMS";flush $TMDout;
   };  
# do we need to get force constant from an existing log file?
   if { $TMDk == -999 } {
    set TMDoutOld [open $TMDlogFileOld r];
    while { [gets $TMDoutOld line] > -1 } {
      set lastline [split $line "\t"];
     }
     close $TMDoutOld
     set TMDk [lindex $lastline 1];
     puts $TMDout "TMD-TCL: Obtained force constant value from file $TMDlogFileOld: $TMDk";flush $TMDout;
   };  
    puts $TMDout "TS \t k \t Force \t Current RMSD \t Target RMSD \t Final RMSD"; flush $TMDout;
}
# now open file for forces output
if { $TMDforcesOutFreq > 0 } {
    if { [file exists $TMDforcesOutFile] } { file copy -force -- "$TMDforcesOutFile" "${TMDforcesOutFile}.BAK"; }
    set TMDforcesOut [open $TMDforcesOutFile w+];
    puts $TMDforcesOut "TMD-TCL: logging forces on $TMDforcedAtoms atoms every $TMDforcesOutFreq steps";
    puts $TMDforcesOut "$TMDforcedAtoms"
    puts $TMDforcesOut "$TMDatoms" ;# list all tmd atoms
    puts $TMDforcesOut "$TMDorientFlags" ;# list all orient atoms
    puts $TMDforcesOut "$TMDforcedFlags" ;# list all forced atoms
    }
###########################################################################################################
proc calcforces { } {

 global TMDatoms TMDorient TMDorientFreq TMDBMD TMDdRMS TMDfirstTstep TMDtargetCoor TMDweights;
 global TMDtargetRMS TMDfinalRMS TMDk TMDfinalK TMDdK TMDout TMDon TMDprintFreq TMDdt TMDdistanceLimit TMDforcedFlags; 
 global TMDorientFlags targetcoor TMDforcesOn;
 global TMDforcesOutFreq TMDforcesOut


 if { $TMDon == 1 } { ;
 
  set currentcoor {};
  loadcoords coor;
  foreach {atom} $TMDatoms { 
   foreach {x y z} $coor($atom) {break;}
   lappend currentcoor $x $y $z ;
  }
 
  set forcelist {}; # used only if forces are output

  set timestep [getstep];
  set dtimestep [expr {$timestep-$TMDfirstTstep}];
  if { ($TMDorient == 1) && ($dtimestep % $TMDorientFreq == 0) } {; # then orient
   set targetcoor [RMSBestFit $TMDtargetCoor $currentcoor $TMDorientFlags]; # use orientation atoms only
  }

  set currentRMS [rmsd $currentcoor $targetcoor $TMDforcedFlags]; # use forced atoms in the computation of force
  if { $TMDtargetRMS == -1 } { set TMDtargetRMS $currentRMS };

  if { $TMDBMD == 0 } {
   if { $TMDtargetRMS > $TMDfinalRMS } { 
    set TMDtargetRMS [expr {$TMDtargetRMS + $TMDdRMS*$TMDdt/1000.0}] 
   } else { 
    set TMDtargetRMS $TMDfinalRMS ; #this will take care of out-of-range values of TargetRMS
   }
  } else {
   if { $currentRMS < $TMDtargetRMS } { set TMDtargetRMS [expr {$TMDtargetRMS + $TMDdRMS}] }; #ratchet down
  }  
  
  if { (( $TMDdK > 0.0 ) &&  ( $TMDk < $TMDfinalK )) || (( $TMDdK < 0.0 ) &&  ( $TMDk > $TMDfinalK )) } {
   set TMDk [expr {$TMDk + $TMDdK*$TMDdt/1000.0}] ; #for gradually increasing/decreasing force constant
  } elseif {$TMDdK !=0.0} {
   set TMDk $TMDfinalK ; # this will take care of out-of-range values of k
  } 
 
  if { ( $TMDforcesOn == 1) && ( $currentRMS > $TMDtargetRMS ) } {; # apply force only if currentRMS above targetRMS
   set forceLog [expr {$TMDk * ($currentRMS - $TMDtargetRMS)}]; # this part of the force due to prefactor only (output)
   set force    [expr {$forceLog/$currentRMS}];                 # this includes normalization

   if { $TMDdistanceLimit > 0.0 } {; #limit maximum force

    if { ( $TMDforcesOutFreq > 0) && ( $timestep % $TMDforcesOutFreq == 0 ) } {; #output forces on TMD atoms

     foreach atom $TMDatoms {xcur ycur zcur} $currentcoor {xtar ytar ztar} $targetcoor forceFlag $TMDforcedFlags {
      set dx [expr {$xtar-$xcur}]; set dy [expr {$ytar-$ycur}]; set dz [expr {$ztar-$zcur}]; 
      set dr [expr {sqrt($dx*$dx+$dy*$dy+$dz*$dz)}];
      set scale [expr { [set scale [expr {$TMDdistanceLimit/$dr} ]] > 1.0 ? 1.0 : $scale }]; # minimum function
      set pre [expr {$force*$forceFlag*$scale}];
      #if dr>TMDdistanceLim, then scale force; i.e this will cap the force! 
      set forcevec [concat [expr {$pre*$dx}] [expr {$pre*$dy}] [expr {$pre*$dz}]];
      addforce $atom $forcevec 
#      lappend forcelist $forcevec
      lappend forcelist [concat $scale $scale $scale];
     };# foreach
    } else {;#dtimestep
     foreach atom $TMDatoms {xcur ycur zcur} $currentcoor {xtar ytar ztar} $targetcoor forceFlag $TMDforcedFlags {
      set dx [expr {$xtar-$xcur}]; set dy [expr {$ytar-$ycur}]; set dz [expr {$ztar-$zcur}]; 
      set dr [expr {sqrt($dx*$dx+$dy*$dy+$dz*$dz)}];
      set scale [expr { [set scale [expr {$TMDdistanceLimit/$dr}]] > 1.0 ? 1.0 : $scale }]; # minimum function
      set pre [expr {$force*$forceFlag*$scale}];
      #if dr>TMDdistanceLim, then scale force; i.e this will cap the force! 
      addforce $atom [concat [expr {$pre*$dx}] [expr {$pre*$dy}] [expr {$pre*$dz}]];
     } ;#foreach
    };#dtimestep
   } else {; #TMDdistanceLimit

    if {( $TMDforcesOutFreq > 0) && ( $timestep % $TMDforcesOutFreq == 0 ) } {; #output forces on TMD atoms

     foreach atom $TMDatoms {xcur ycur zcur} $currentcoor {xtar ytar ztar} $targetcoor forceFlag $TMDforcedFlags {
      set pre [expr {$force*$forceFlag}];
      set forcevec [concat [expr {$pre*($xtar-$xcur)}] [expr {$pre*($ytar-$ycur)}] [expr {$pre*($ztar-$zcur)}]];
      addforce $atom $forcevec
      lappend forcelist $forcevec
     };#foreach
    } else {;#dtimestep
     foreach atom $TMDatoms {xcur ycur zcur} $currentcoor {xtar ytar ztar} $targetcoor forceFlag $TMDforcedFlags {
      set pre [expr {$force*$forceFlag}];
      addforce $atom [concat [expr {$pre*($xtar-$xcur)}] [expr {$pre*($ytar-$ycur)}] [expr {$pre*($ztar-$zcur)}]];
     };#foreach
    }; #dtimestep
    
   };#if TMDdistanceLimit
#  output forces
   if { ( $TMDforcesOutFreq > 0) && ( $timestep % $TMDforcesOutFreq == 0 ) } { ; #output forces
    puts $TMDforcesOut "$timestep \t $forcelist"; flush $TMDforcesOut;
   }

  } else {; #TMDforcesOn
    set forceLog 0.0
  }  ; #apply force
  if { $dtimestep % $TMDprintFreq == 0 } { ; #print
   puts $TMDout "$timestep \t $TMDk \t $forceLog \t $currentRMS \t $TMDtargetRMS \t $TMDfinalRMS"; flush $TMDout;
  }
  
 } ;# TMDon 
} ;#calcforces
