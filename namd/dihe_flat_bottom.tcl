# adopted from NAMD 2.7b user guide
# The IDs of the four atoms defining the phi-dihedral
# see Ovchinnikov & Karplus 2013 (SCM paper)
# the potential is U = Kphi * max ( 0 , |phi(X) - phi_0| - Delta_phi )^2 / 2
  set PI 3.14159265359879 ;
#
  set aid1 [atomid diala 1 C];
  set aid2 [atomid diala 2 N];
  set aid3 [atomid diala 2 CA];
  set aid4 [atomid diala 2 C];
#  set phi0 -115 ;# c7eq
#  set dphi0 115 ;# c7eq
  set phi0 65 ;# c7ax
  set dphi0 65 ;# c7ax
#
  set phi0  [expr {  $phi0 * $PI / 180 } ]; # convert to radians
  set dphi0 [expr { $dphi0 * $PI / 180 } ];
#
# Spring constant
  set Kphi 20.0 ; # kcal/mol/rad^2 
#
  addatom $aid1
  addatom $aid2
  addatom $aid3
  addatom $aid4
#
  proc calcforces {} {
    global aid1 aid2 aid3 aid4 Kphi phi0 dphi0 PI
    loadcoords r
# Calculate the current dihedral
    set phi [getdihedral $r($aid1) $r($aid2) $r($aid3) $r($aid4) ]; # in degrees
    set phi [expr {$phi * $PI / 180.0 } ] ; # radians
# ensure dihedral within [-180 ... 180 ]
    if { $phi > $PI } { set phi [ expr { $phi - 2.0*$PI } ] } elseif { $phi < -$PI } { set phi [ expr { $phi + 2.0*$PI } ] }
# dihedral difference from center of basin:
    set dphi [expr { $phi - $phi0 }];
# ensure (again) dihedral within [-180 ... 180 ]
    if { $dphi > $PI } { set dphi [ expr { $dphi - 2.0*$PI } ] } elseif { $dphi < -$PI } { set dphi [ expr { $dphi + 2.0*$PI } ] }
# check if dihedral outside of allowed region
    set d [expr { abs($dphi) - $dphi0 }];
    if { $d > 0 } {
# (optional) Add restraining energy to MISC in the energy output
     addenergy [expr { 0.5 * $Kphi * $d * $d } ];
# Calculate the force on the dihedral according to the harmonic restraint
     set force [expr { -$Kphi * $d } ];
     if { $dphi < 0 } { set force [expr { -$force } ] }
# Calculate the gradients
     foreach {g1 g2 g3 g4} [dihedralgrad $r($aid1) $r($aid2) $r($aid3) $r($aid4) ] {}
# The force to be applied on each atom is proportional to its
# corresponding gradient
     addforce $aid1 [vecscale $g1 $force ]
     addforce $aid2 [vecscale $g2 $force ]
     addforce $aid3 [vecscale $g3 $force ]
     addforce $aid4 [vecscale $g4 $force ]
   } ;# d>0
 } ;# calcforces
