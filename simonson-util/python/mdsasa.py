#!/bin/python2
from __future__ import print_function
from sys import stdout, stderr, exit, argv
import mdtraj as md

if len(argv) < 3 :
 raise ValueError("USAGE : ./mdsasa.py <psf_file> <dcd_file>")

myname=argv[0];
psffile=argv[1];
dcdfile=argv[2];
if len(argv) > 3 :
 stride=int(argv[3]);
#
# ==============================================
# ========================== Aux print functions
def dprint(*args):
  print(" ===>",myname,": ",end="");
  for arg in args:
   print(arg,end="")
  print(""); # flush
#==========================
def derror(*args):
  print(" ===>",myname,"ERROR : ",end="");
  for arg in args:
   print(arg,end="")
  print(""); # flush
# ==========================
dprint("Reading simulation coordinates from file '",dcdfile,"'");
try :
 stride
except NameError:
 stride=1 ;# read all by default
#
dcd=md.load_dcd(dcdfile,top=psffile, stride=stride);
dprint("Computing SASA for ", dcd.n_frames," frames:" );
sasa=md.shrake_rupley(dcd, mode="residue");
total_sasa=sasa.sum(axis=1); # sum over residues
for s in total_sasa :
 print("%12.5f" %s);
