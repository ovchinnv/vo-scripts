source userani.vmd

add_section kickoff -start first
add_section showstart -start first -end 5% -noplay
add_section anirasterfoc -after showstart -shiftend "7%"
add_section aniraster -like anirasterfoc -shiftstart "1%"
add_section staticraster -after aniraster -end last

add_eval kickoff display resetview
add_eval kickoff material change opacity Glossytrans1 1.000000
add_eval kickoff rotate x by 180
add_eval showstart dummy
add_userani aniraster raster_replicas 40 40
add_userani anirasterfoc AutoFocusAllVisible
add_userani anirasterfoc AutoScaleAllVisible
add_eval staticraster raster_replicas 40 40

#more complex example
#fly and zoom from replica to replica among the raster
set each 4%
set each2 [expr% int($each/2.)]
for {set i 0} {$i < $num_replicas} {incr i} {
  if {$i == 0} {
    add_section focus${i} -after aniraster -shiftend $each
    add_section autoscale${i} -after aniraster -shiftend $each2
    add_section zoomout${i} -like autoscale${i} -shiftall $each2
  } else {
    add_section focus${i} -like focus[expr $i-1] -shiftall $each
    add_section autoscale${i} -like autoscale[expr $i-1] -shiftall $each
    add_section zoomout${i} -like zoomout[expr $i-1] -shiftall $each 
  }

  if {$i < [expr $num_replicas -1]} {
    add_userani zoomout${i} FadeZoomFromCurrent 0.5
  }
  add_userani focus${i} AutoFocusMolidsVisible $i
  add_userani autoscale${i} AutoScaleMolidsVisible $i
}

add_section restauto -after "focus[expr $num_replicas-1]" -shiftend 10%
add_userani restauto AutoFocusAllVisible
add_userani restauto AutoScaleAllVisible

add_section focus10fadeothers -after restauto -shiftstart -5% -shiftend 5%
add_section focus10 -after focus10fadeothers -shiftend 3%

add_userani focus10 AutoFocusMolidsVisible 10
add_userani focus10 AutoScaleMolidsVisible 10
add_userani focus10fadeothers FadeTransparency Glossytrans1 1 0

#just print what we have setup on source
print_story
print_timeline 2

