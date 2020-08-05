#!/bin/vmd
set name [exec ls | grep _msn.pdb | awk -F "_" "{print \$1}"]
puts $name
