#!/bin/vmd
set name [exec ls | grep _msn.pdb | awk -F "_msn" "{print \$1}"]
puts $name
