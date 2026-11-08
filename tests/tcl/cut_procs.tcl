# Cuts procedures out of the original's AppMain.tcl, unchanged, for the drivers in this folder.
# usage: tclsh cut_procs.tcl /Applications/ApE.app/Contents/Resources/Scripts/AppMain.tcl name name ... > orig.tcl
# A procedure is cut from its `proc name ` line (the last one in the file, as that is the one that takes effect when the file
# loads) to the line where the text is a complete command (`info complete`), so no line numbers have to be kept.
# The original's own code is not part of this repository: cut it from your copy of the program.
set f [open [lindex $argv 0] r]
set lines [split [read $f] "\n"]
close $f
foreach name [lrange $argv 1 end] {
  set start -1
  for {set i 0} {$i < [llength $lines]} {incr i} {
    if {[string first "proc $name " [lindex $lines $i]] == 0} { set start $i }
  }
  if {$start < 0} { puts stderr "no proc $name"; exit 1 }
  set text ""
  for {set j $start} {$j < [llength $lines]} {incr j} {
    append text [lindex $lines $j] "\n"
    if {[info complete $text]} { break }
  }
  puts $text
}
