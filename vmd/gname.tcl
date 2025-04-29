#!/bin/vmd
#set name [exec ls | grep -E sn.pdb | sed -E s/_m?sn.pdb//]
set name [exec ls | grep -E sn.psf | sed -E s/_sn.psf//]
puts $name
