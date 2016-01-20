proc makemov { {fname "movie"} {id 0} {renderer "snapshot"} {freq 1} } {
        #{ mama {id 0} {renderer Tachyon} {freq 1} {fname movie\/movie} }
	# get the number of frames in the movie
	set num [molinfo $id get numframes]
	# loop through the frames
	for {set i 0} {$i < $num} {incr i $freq} {
		animate goto $i
                display update
		set filename $fname.[format "%04d" [expr $i/$freq]]
		if { $renderer == "Tachyon" } {
		  render $renderer $filename.dat
		  exec /usr/local/bin/tachyon -format RGB $filename.dat -o $filename.rgb  
		  exec convert -quality 100 $filename.rgb $filename.jpg
                  exec rm -f $filename.dat
		  exec rm $filename.rgb
                }
		if {$renderer == "snapshot"} {
		   render $renderer $filename.rgb
                   exec convert -quality 100 $filename.rgb $filename.jpg
                   exec rm -f $filename.dat
		   exec rm $filename.rgb
		 } 
         }	         
}
#makemov "movie\/myorp2" [molinfo top]

#makemov "movie2/myo"


