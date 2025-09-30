#!/bin/python
from pymol import pm

pm.delete('all')
pm.set('stereo','on')
pm.stereo('walleye')

#names=[ 'h7vrcc15st_A', 'fluab_HA1', 'fluab_HA2' ]
names=[ 'h7vrcc15st_A', 'fluab_HA' ]
basename='aligned';
#names=[ 'h7vrcc15st_A' ];
colors=[ 'red', 'green', 'blue'];
i=0;
for name in names:
 obj='aligned'+name; # these will become objects for pymol selections
 pdb=obj+'.pdb';
 print(pdb)
 pm.load(pdb)
# pm.hide(obj)
 pm.hide('lines',obj) ;# should work similarly for all selections
# pm.show('spheres',obj) ;# this works
 pm.show('cartoon',obj) ;# this works
# pm.show('ribbon',obj) ;# ribbon rep is crap
# pm.set('cartoon_color',colors[i]); # do not use like thism because it will set all cartoons to the same color
 col=colors[i];
 pm.color(col,obj)
 i=i+1

pm.set('cartoon_fancy_helices')
pm.set('cartoon_highlight_color','gray75')

# alignments
#aln1=pm.align(basename+'h7vrcc15st_A', basename+'fluab_HA1') ;#
aln1=pm.align(basename+'h7vrcc15st_A', basename+'fluab_HA') ;#
#aln1=pm.align(basename+'h7vrcc15st_A////CA', basename+'fluab_HA1////CA') ;# this is a lower RMSD
print('RMSD:' + str(aln1[0])+' based on ' + str(aln1[1]) + ' atoms')
#print(f"{aln1[0]}"); # incompatible, don't know why
#aln2=pm.super(basename+'h7vrcc15st_A', basename+'fluab_HA1') ;# slightly lower RMSD than aln1, fewer atoms
#print("RMSD:" + str(aln2[0])+" based on " + str(aln2[1]) + " atoms") ;#
#
#aln2=pm.align(basename+'h7vrcc15st_A', basename+'fluab_HA2') ;#
#print("RMSD:" + str(aln2[0])+" based on " + str(aln2[1]) + " atoms")

