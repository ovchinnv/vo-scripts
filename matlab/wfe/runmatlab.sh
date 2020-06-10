#!/bin/bash
# Simone's driver script to run wfet.m on all mutants
#
#module load matlab

#cd titan

for dir in `ls -d */`; do

    echo $dir

    cd $dir

    for win in `seq 1 7`; do
        echo -n > fbwin_${win}.dat
#        cp fbwin2_${win}.dat fbwin_${win}.dat
        for step in `seq 1 2`; do
            cat fbwin${step}_${win}.dat >> fbwin_${win}.dat
        done
        wc -l fbwin_${win}.dat
    done


    #cat fbwin2_1.dat.dat fbwin3_1.dat fbwin4_1.dat fbwin5_1.dat > fbwin_1.dat 
    #cat fbwin2_2.dat.dat fbwin3_2.dat fbwin4_2.dat fbwin5_2.dat > fbwin_2.dat
    #cat fbwin2_3.dat.dat fbwin3_3.dat fbwin4_3.dat fbwin5_3.dat > fbwin_3.dat
    #cat fbwin2_4.dat.dat fbwin3_4.dat fbwin4_4.dat fbwin5_4.dat > fbwin_4.dat
    #cat fbwin2_5.dat.dat fbwin3_5.dat fbwin4_5.dat fbwin5_5.dat > fbwin_5.dat
    #cat fbwin2_6.dat.dat fbwin3_6.dat fbwin4_6.dat fbwin5_6.dat > fbwin_6.dat
    #cat fbwin2_7.dat.dat fbwin3_7.dat fbwin4_7.dat fbwin5_7.dat > fbwin_7.dat

    ln -sf ../wfet.m 
    #matlab -nodisplay -nodesktop -r "run wfet.m"
    octave < wfet.m
    cp dg.txt ../dg_${dir%/}.dat
    cp wfet.eps ../fet_${dir%/}.eps
    cd ../
#    exit
done

