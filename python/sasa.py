#!/bin/python
from __future__ import print_function
from sys import stdout, stderr, exit, argv
import mdtraj as md

print("%",md.__file__)

myname=argv[0];
# ========================== Aux print functions
def dprint(*args):
  print("% ===>",myname,": ",end="");
  for arg in args:
   print(arg,end="")
  print(""); # flush
#==========================
def derror(*args):
  print("% ===>",myname,"ERROR : ",end="");
  for arg in args:
   print(arg,end="")
  print(""); # flush
# ==========================

if len(argv) < 2 :
 raise ValueError("USAGE : ./mdsasa.py <struc_file> <dcd_file>(optional)")
else:
 topfile=argv[1];
 if len(argv) < 3 :
  qdcd=0
 else:
  qdcd=1;
  dcdfile=argv[2];
  if len(argv) > 3:
   stride=int(argv[3]); # only makes sense if dcd present
  else:
   stride=1;
# ==============================================
struc=md.load(topfile);

# ==============================================
if (qdcd):
 dprint("Reading simulation coordinates from file '",dcdfile,"'");
 dcd=md.load_dcd(dcdfile,stride=stride);
#
nm2a=100 ; #convert to A from nm
if struc.n_frames>0 :
 dprint("Computing SASA for ", struc.n_frames," frames:" );
 sasa=md.shrake_rupley(struc, mode="residue"); # mode residue means one number per resdiue (rather than mode="atom"

for i in range(sasa.shape[0]):
 dprint("Frame ",i);
 for s in sasa[i]:
  print("%12.6f" %(s*nm2a));

 dprint("Total area")
 total_sasa=sasa[i].sum(); # sum over residues
 print('%%%12.5f' %(total_sasa*nm2a));

