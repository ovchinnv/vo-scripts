#!/bin/vmd
set name [exec ls | grep -E sn.pdb | sed -E s/_m?sn.pdb//]
puts $name
