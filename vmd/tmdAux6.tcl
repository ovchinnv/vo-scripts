#!/bin/tcl
# version 3: improved efficiency
# version 5: bugfixes: (1) eigenvalues were not sorted (2) reflections were allowed in best fit
# version 6: bugfixes
# Victor Ovchinnikov, Harvard University/MIT, 2007
################################## subroutines ##################
proc multMat {m1 m2} {
 set result {};
 foreach {b11 b12 b13 b21 b22 b23 b31 b32 b33} $m2 {break};
 foreach {a1 a2 a3} $m1 { 
  lappend result [expr {$a1*$b11+$a2*$b21+$a3*$b31}] [expr {$a1*$b12+$a2*$b22+$a3*$b32}] [expr {$a1*$b13+$a2*$b23+$a3*$b33}];
 };
 return $result 
}

proc multTrMat {m1 m2} { ;# does m1 x tr(m2)
 set result {};
 foreach {a1 a2 a3} $m1 {
  foreach {b1 b2 b3} $m2 {
   lappend result [expr {$a1*$b1+$a2*$b2+$a3*$b3}]
  }
 }
 return $result
}
proc multTrMat2 {m1 m2} { ;# does tr(m1) x m2
 set m1 [transpose $m1];
 set m2 [transpose $m2];
 return [multTrMat $m1 $m2];
}

proc transpose {m1} {
 foreach {b11 b12 b13 b21 b22 b23 b31 b32 b33} $m1 {break};
 return [concat $b11 $b21 $b31 $b12 $b22 $b32 $b13 $b23 $b33];
} 

proc vecdot {v1 v2} {
 foreach {a1 a2 a3} $v1 {b1 b2 b3} $v2 {break};
 return [expr {$a1*$b1+$a2*$b2+$a3*$b3}];
}  
proc veccross {v1 v2} {
 foreach {a1 a2 a3} $v1 {b1 b2 b3} $v2 {break};
 return [concat [expr {$a2*$b3-$a3*$b2}] [expr {$a3*$b1-$a1*$b3}] [expr {$a1*$b2-$a2*$b1}]];
}  
proc unit {v1} { ;# normalize vector to unity
 set errTol 1.0e-9
 foreach {x y z} $v1 {break;}
  set n [expr {sqrt($x*$x+$y*$y+$z*$z)} ];
  if { [expr {$n<$errTol}] } {set n 1.0} ;# protect against divide overflow
  return [concat [expr {$x/$n}] [expr {$y/$n}] [expr {$z/$n}] ];
}

proc vecnorm {v1} { ;# compute norm of vector
 foreach {x y z} $v1 {break;}
  return [expr {sqrt($x*$x+$y*$y+$z*$z)} ];
}

proc vecsub {v1 v2} { 
 foreach {a1 a2 a3} $v1 {b1 b2 b3} $v2 {break};
 return [concat [expr {$a1-$b1}] [expr {$a2-$b2}] [expr {$a3-$b3}] ];
}

proc vecscale {s v1} { ;
 foreach {a1 a2 a3} $v1 {break};
 return [ concat [expr {$s*$a1}] [expr {$s*$a2}] [expr {$s*$a3}] ];
}

proc rmsd {v1 v2 weights} {
 set rmsd 0.0;
 set natom 0.0;
 foreach {a b c} $v1 {d e f} $v2 w $weights {
# this line (rmsd=...)fails on daedalus1:
#  set dx [expr {$a-$d}]
#  set dy [expr {$b-$e}]
#  set dz [expr {$c-$f}]
  set rmsd [expr {$rmsd+$w*(($a-$d)*($a-$d) + ($b-$e)*($b-$e) + ($c-$f)*($c-$f))}];
#  set rmsd [expr {$rmsd + $w * ( $dx*$dx + $dy*$dy + $dz*$dz ) } ];
  set natom [expr { $natom + $w }];
  }
  return [expr {sqrt($rmsd/$natom)}];
}

proc moveToO {v1 masses} { ; # move v1 to origin & compute COM  
 set xCOM 0.0; set yCOM 0.0; set zCOM 0.0;
 set totmass 0.0
 foreach {x1 y1 z1} $v1 mass $masses { 
  set xCOM [expr {$xCOM + $mass*$x1}];
  set yCOM [expr {$yCOM + $mass*$y1}];
  set zCOM [expr {$zCOM + $mass*$z1}];
  set totmass [expr {$totmass+$mass}];
 }
 set xCOM [expr {$xCOM/$totmass}];
 set yCOM [expr {$yCOM/$totmass}];
 set zCOM [expr {$zCOM/$totmass}];
 set v3 {};
 foreach {x1 y1 z1} $v1 { lappend v3 [expr {$x1-$xCOM}] [expr {$y1-$yCOM}] [expr {$z1-$zCOM}] }; # translated coords.
 return [concat $v3 ":" $xCOM $yCOM $zCOM]
} 

proc translate {d v1} { ; # move v1 to origin & compute COM  
 foreach {dx dy dz} $d {break}; 
 set v3 {};
 foreach {x1 y1 z1} $v1 { 
  lappend v3 [expr {$x1+$dx}] [expr {$y1+$dy}] [expr {$z1+$dz}] 
 } 
 return $v3
} 

proc det {m} {
 foreach {a11 a12 a13 a21 a22 a23 a31 a32 a33} $m {break;}
 return [expr $a11*$a22*$a33+$a12*$a23*$a31+$a13*$a21*$a32-$a13*$a22*$a31-$a12*$a21*$a33-$a11*$a23*$a32];
}
##################################################################
proc eig {M} { ;# Exact diagonalization of a 3X3 SYMMETRIC matrix

set pi 3.14159265358979;
set errTol 1.0e-09;
#puts "$M"
# get components
foreach {a11 a12 a13 a21 a22 a23 a31 a32 a33} $M {break}; 
# solve the characteristic equation
set a2 [expr {-($a11+$a22+$a33)}];
set a1 [expr { ($a11*$a22+$a11*$a33+$a22*$a33-$a23*$a32-$a12*$a21-$a13*$a31)}];
set a0 [expr {-($a11*$a22*$a33+$a12*$a23*$a31+$a13*$a21*$a32-$a11*$a23*$a32-$a13*$a31*$a22-$a12*$a21*$a33)}];

set Q [expr {(3.0*$a1-$a2*$a2)/9.0}];
set R [expr {(9.0*$a2*$a1-27.0*$a0-2*$a2*$a2*$a2)/54.0}];
set D [expr {$Q*$Q*$Q+$R*$R}];

set rootQ [expr {sqrt(-$Q)}];
set rootQ3 [expr {$rootQ*$rootQ*$rootQ}];

if { $rootQ == 0 } { ;# for the indentity matrix
 set z1 [expr {-$a2/3.0}];
 set z2 [expr {-$a2/3.0}];
 set z3 [expr {-$a2/3.0}];
} else {
 set theta [expr {acos($R/$rootQ3)}];
 set z1 [expr {2.0*$rootQ*cos($theta/3.0)-$a2/3.0}];
 set z2 [expr {2.0*$rootQ*cos(($theta+2.*$pi)/3.)-$a2/3.0}];
 set z3 [expr {2.0*$rootQ*cos(($theta+4.*$pi)/3.)-$a2/3.}];
}
############# truncate if within errTol of 0 ############################
if { [expr {abs($z1) < $errTol }] } { set z1 0. };
if { [expr {abs($z2) < $errTol }] } { set z2 0. };
if { [expr {abs($z3) < $errTol }] } { set z3 0. };
set ev [lsort -real [concat $z1 $z2 $z3] ];
######################### now compute eigenvectors #################################
set evec {};

foreach mu $ev {
 foreach {a b c d e f g h i} [concat [expr {$a11-$mu}] $a12 $a13 $a21 [expr {$a22-$mu}] $a23 $a31 $a32 [expr {$a33-$mu}]] {break;} ;# notation change
 set P [unit [concat $a $d $g]];
 set Q [unit [concat $b $e $h]];
 set R [unit [concat $c $f $i]];

#puts "$a $b $c"
#puts "$d $e $f"
#puts "$g $h $i"
#puts "$mu $P [expr {abs([vecnorm $P])}]"

 if { [expr {abs([vecnorm $P])}] < $errTol } { set v {1.0 0.0 0.0}
 } elseif { [expr {abs([vecnorm $Q])}] < $errTol } { set v {0.0 1.0 0.0}
 } elseif { [expr {abs([vecnorm $R])}] < $errTol } { set v {0.0 0.0 1.0}
 } elseif { [expr {abs(abs([vecdot $P $Q]) - 1.0)  > $errTol}] } {; # P || Q
 # z=1
  set det12 [expr {$a*$e-$b*$d}];
  set det13 [expr {$a*$h-$b*$g}];

  if { [expr {abs($det12) > $errTol}] } {
   set v [unit [concat [expr {($b*$f-$c*$e)/$det12}] [expr {($d*$c-$a*$f)/$det12}] 1.0]];
  } elseif { [expr {abs($det13) >$errTol }]} { 
   set v [unit [concat [expr {($b*$i-$c*$h)/$det13}] [expr {($g*$c-$a*$i)/$det13}] 1.0]];
  }
 } else {
   # otherwise there is a vector with z = 0; try setting x=1
  if { [expr {abs($b) > $errTol}] } { 
    set v [unit [concat 1.0 [expr {-$a/$b}] 0.0 ]]; 
  } elseif { [expr {abs($e) > $errTol}] } { 
    set v [unit [concat 1.0 [expr {-$d/$e}] 0.0 ]]; 
  } elseif { [expr {abs($h) > $errTol} ] } { 
    set v [unit [concat 1.0 [expr {-$g/$h}] 0.0]];
  } else { ;# Q=0; this should not be possible because P Q R are checked above !
    set v {0.0 1.0 0.0}
  }
 } 
# puts "$mu $v"
 lappend evec $v;
}  

#DONE : output
foreach {v1 v2 v3} $evec {break;}
# make sure that triple degeneracy is taken into account
if {[expr {abs(abs([vecdot $v1 $v2]) - 1.0) < $errTol}] & [expr {abs(abs([vecdot $v2 $v3]) - 1.0) < $errTol}]} {;#
 set v1 {1.0 0.0 0.0};
 set v2 {0.0 1.0 0.0};
 set v3 {0.0 0.0 1.0};
} 
# make sure that double degeneracy is taken into account
if {[expr {abs(abs([vecdot $v1 $v2]) - 1.0) < $errTol}] } {set v2 [veccross $v3 $v1]; puts $v2}
if {[expr {abs(abs([vecdot $v2 $v3]) - 1.0) < $errTol}] } {set v3 [veccross $v1 $v2]; puts $v3}
if {[expr {abs(abs([vecdot $v3 $v1]) - 1.0) < $errTol}] } {set v1 [veccross $v2 $v3]; puts $v1}

set evec [transpose [ concat $v1 $v2 $v3 ]]; # this is the matrix of eigen-vectors
#puts $evec
#puts "[multTrMat2 $evec $evec]" ;# this should give the identity matrix
#puts "EIGEN: $ev"
#foreach {x y z} $evec {puts "EIGEN $x $y $z"};
return [concat $evec ":" $ev] ;# append eigenvalues at the end
}
###############################################################################################
proc RMSBestFit {x0 y0 weights} { ;# implementation of the algorithm by Kabsch 1976
foreach {x} [split [moveToO $x0 $weights] : ] {break} ; # compute translate to origin
foreach {y COMy} [split [moveToO $y0 $weights] : ] {break} ; # compute COM & translate to origin

 foreach {r11 r12 r13 r21 r22 r23 r31 r32 r33} {0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0} {break}; # initialize
 foreach {x1 x2 x3} $x {y1 y2 y3} $y {w} $weights {
  set r11 [expr {$r11+$y1*$x1*$w}];
  set r12 [expr {$r12+$y1*$x2*$w}];
  set r13 [expr {$r13+$y1*$x3*$w}];
#
  set r21 [expr {$r21+$y2*$x1*$w}];
  set r22 [expr {$r22+$y2*$x2*$w}];
  set r23 [expr {$r23+$y2*$x3*$w}];
#
  set r31 [expr {$r31+$y3*$x1*$w}];
  set r32 [expr {$r32+$y3*$x2*$w}];
  set r33 [expr {$r33+$y3*$x3*$w}];
 }
 set r [concat $r11 $r12 $r13 $r21 $r22 $r23 $r31 $r32 $r33];
 set rr [multTrMat2 $r $r]; # r'r
 set detr [det $r]; # determinant of r; if this is negative, we will have a reflection;
 set reflect_flag [expr {$detr/abs($detr)}]; # this flag is -1 if we have a reflection
# let's diagonalize rr
 foreach {a mu} [split [eig $rr] : ] {break}; # separate evectors & evalues
 set b [multMat $r $a];

 foreach {mu1 mu2 mu3} $mu {break}; 
 if { $mu1 == 0.0 } {set mu1 1.0}; if { $mu2 == 0.0 } {set mu2 1.0}; if { $mu3 == 0.0 } {set mu3 1.0};
 set MUinv  [concat [expr {$reflect_flag/sqrt($mu1)}] 0. 0. 0. [expr {1.0/sqrt($mu2)}] 0. 0. 0. [expr {1.0/sqrt($mu3)}]]; # eigenvalue matrix
 set b [multTrMat $b $MUinv];
 set u [multTrMat $b $a]; # this is the rotation matrix
 set x [translate $COMy [multTrMat $x $u]]; # apply rotation & move to the center of mass of y
# puts "\n";foreach {q v w} $mu {puts "$q $v $w"}
 return $x ; # return the new vector
} 
