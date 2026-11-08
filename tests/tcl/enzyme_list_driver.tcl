# Runs the original's List Enzymes and Quick Lists code, unchanged (orig.tcl is cut out of AppMain.tcl with cut_procs.tcl),
# against a stand-in text window, and prints the cases and what the original wrote into its result window.
# This is how src/seq/enzlist/golden.txt was made (see MECHANICS.md, section 8, "List Enzymes and Quick Lists").
#
# usage: tclsh cut_procs.tcl AppMain.tcl revcom validation_code read_enzymes not_in_list initialize_methylation \
#          toggle_methylation findenzymes enz_site_usage_filter enz_count list_enzymes_window list_enzymes_text > orig.tcl
#        tclsh enzyme_list_driver.tcl orig.tcl <seed> <cases> > golden.txt          (committed: seed 20261107, 72 cases)
#        tclsh enzyme_list_driver.tcl orig.tcl <seed> <cases> 1 > golden.txt        (with runs that make sites overlap)
# A larger file than the committed one can be checked with the ignored test in src/seq/enzlist/tests.rs.
# The original's enzyme file (Default_Enzymes.txt) is read by the original's own `read_enzymes`; its sites are found by the
# original's own `findenzymes`, with Dam/Dcm methylation on or off as the case says.
set origfile [lindex $argv 0]
set seed [lindex $argv 1]
set ncases [lindex $argv 2]
# with a fourth argument of 1, runs of a site repeated at its own period are put in, so that sites overlap
set dense [expr {[lindex $argv 3] eq "1"}]
set ACC "/Applications/ApE.app/Contents/Resources/Scripts/Accessory Files"

# ---- stand-ins for the windowing code
proc mc {args} { if {[llength $args] > 1} { return [format [lindex $args 0] {*}[lrange $args 1 end]] }; return [lindex $args 0] }
proc tk_messageBox {args} { return ok }
proc warning_dialog {args} {}
proc add_compatibles {} {}
proc wm {args} { if {[lindex $args 0] eq "title"} { return TITLE }; return "" }
proc bind {args} {}
proc bindtags {args} {}
proc window_on_screen {args} {}
proc focus {args} {}
proc textbox_geometry {args} { return "" }
proc winfo {args} { return 0 }
proc add_output_menu {args} { return "" }
proc create_textframe {a h} {}
proc new_analysis_window {w title name} { return .an1 }
proc bp2ix {t n} { return $n }
proc end_ix {t} { return end }
proc ix2bp {t ix} { if {$ix eq "end"} { return [string length $::seq] }; if {$ix eq "1.0"} { return 0 }; return $ix }
proc textarea_get {t i1 i2} { set a [ix2bp $t $i1]; set b [ix2bp $t $i2]; return [string range $::seq $a [expr {$b-1}]] }
proc .w1.textarea {cmd args} {
  switch -- $cmd {
    tag { if {[lindex $args 0] eq "ranges" && [lindex $args 1] eq "sel"} { return $::sel }; return {} }
    index { if {[lindex $args 0] eq "sel.first"} { return [lindex $::sel 0] }; if {[lindex $args 0] eq "sel.last"} { return [lindex $::sel 1] } }
  }
}
proc .an1.textframe.text {cmd args} {
  if {$cmd ne "insert"} { return }
  append ::out [lindex $args 1]
}

source $origfile

set w .w1
set info(Accdir) $ACC
set info(clock_format) "DATE"
set info(last_analysis_xy) {0 0}
set dialogblock 0
set ok 0
read_enzymes [file join $ACC Enzymes Default_Enzymes.txt]

expr {srand($seed)}
proc rnd {n} { return [expr {int(rand()*$n)}] }
proc randdna {n} { set s ""; for {set i 0} {$i < $n} {incr i} { append s [string index ACGT [rnd 4]] }; return $s }
proc concrete {pat} {
  set out ""
  foreach c [split $pat ""] {
    set bases [string map {A A C C G G T T R AG Y CT S CG W AT K GT M AC B CGT D AGT H ACT V ACG N ACGT} $c]
    append out [string index $bases [rnd [string length $bases]]]
  }
  return $out
}
# a random sequence with the sites of several enzymes put in (on either strand, some twice), so that the lists have
# enzymes with no site, one site and several
proc make_sequence {len} {
  set s [randdna $len]
  set n [expr {4 + [rnd 14]}]
  for {set k 0} {$k < $n} {incr k} {
    set e [lindex $::enzymes [rnd [llength $::enzymes]]]
    set site [concrete $::enzinfo(flatpat,$e)]
    if {[rnd 2]} { set site [revcom $site] }
    if {[string length $site] >= $len} { continue }
    set at [rnd [expr {$len - [string length $site]}]]
    set s [string replace $s $at [expr {$at + [string length $site] - 1}] $site]
    if {[rnd 4] == 0} {
      set at [rnd [expr {$len - [string length $site]}]]
      set s [string replace $s $at [expr {$at + [string length $site] - 1}] $site]
    }
  }
  if {$::dense} {
    for {set k 0} {$k < 2 + [rnd 3]} {incr k} { set s [put_run $s] }
  }
  return $s
}
# a site repeated at its own period (the smallest shift that lays it over itself), so that its copies overlap as far as
# it allows, put in at a random place, on either strand
proc put_run {s} {
  set e [lindex $::enzymes [rnd [llength $::enzymes]]]
  set site [concrete $::enzinfo(flatpat,$e)]
  if {[rnd 2]} { set site [revcom $site] }
  set n [string length $site]
  set period $n
  for {set p 1} {$p < $n} {incr p} {
    if {[string range $site $p end] eq [string range $site 0 [expr {$n - $p - 1}]]} { set period $p; break }
  }
  set copies [expr {2 + [rnd 4]}]
  set total [expr {$n + $period * ($copies - 1)}]
  set run [string range [string repeat [string range $site 0 [expr {$period - 1}]] [expr {$total / $period + 2}]] 0 [expr {$total - 1}]]
  if {[string length $run] >= [string length $s]} { return $s }
  set at [rnd [expr {[string length $s] - [string length $run]}]]
  return [string replace $s $at [expr {$at + [string length $run] - 1}] $run]
}
proc pick_some {n} {
  set chosen [list]
  while {[llength $chosen] < $n} {
    set e [lindex $::enzymes [rnd [llength $::enzymes]]]
    if {$e ni $chosen} { lappend chosen $e }
  }
  return [lsort -unique $chosen]
}

set whiches {all selected absent unique}
for {set k 1} {$k <= $ncases} {incr k} {
  set len [expr {50 + [rnd 190]}]
  set circ [rnd 2]
  set methylated [rnd 2]
  set ::seq [make_sequence $len]
  set ::sel {}
  if {[rnd 3] > 0} {
    set a [rnd [expr {$len - 20}]]
    set b [expr {$a + 10 + [rnd [expr {$len - $a - 9}]]}]
    set ::sel [list $a $b]
  }
  set which [lindex $whiches [expr {$k % 4}]]
  set format [rnd 4]
  set columns [expr {1 + [rnd 10]}]
  set region [expr {[rnd 3] == 0 ? "all" : "selection"}]
  set chosen [expr {$k % 9 == 0 ? {} : [pick_some [rnd 14]]}]
  set info(enz_currently_selected) $chosen
  set info(enz_list_columns) $columns
  set info($w,circular) [expr {$circ ? "circular" : "linear"}]
  set info($w,scanned) 0
  set info($w,Dam_Dcm_methylated) $methylated
  findenzymes $w
  puts "CASE list $k"
  puts "circular $circ"
  puts "seq $::seq"
  puts "sel [expr {$::sel eq {} ? "-" : $::sel}]"
  puts "enz [expr {[llength $chosen] == 0 ? "-" : $chosen}]"
  puts "methylated $methylated"
  puts "format $format"
  puts "which $which"
  puts "columns $columns"
  puts "region $region"
  foreach e $::enzymes {
    if {[llength $info($w,$e)] > 0} { puts "sites $e [join $info($w,$e) { }]" }
  }
  set ::out ""
  list_enzymes_window $w $format $which $region
  puts "out"
  foreach line [split $::out "\n"] { puts "|$line" }
  puts "END"
}
