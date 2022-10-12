#!/bin/vmd
set name [exec ls | grep sn.pdb | sed -E s/_m?sn.pdb//]
puts $name
