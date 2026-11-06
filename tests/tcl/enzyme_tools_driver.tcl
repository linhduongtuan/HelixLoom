# Runs the original's Silent Sites and Find Diagnostic Site code, unchanged (orig.tcl is cut out of AppMain.tcl),
# against a stand-in text window, and prints the cases and what the original wrote into its result window.
# This is how src/seq/silent/golden.txt and src/seq/diagnostic/golden.txt were made (see MECHANICS.md, section 8).
#
# usage: tclsh enzyme_tools_driver.tcl orig.tcl <seed> <silent-cases> <diag-cases> > golden_all.txt
#   (the committed files: seed 20261106, 80 silent cases, 30 diagnostic cases; the part from `CASE silent` up to the first
#   `CASE diagnostic` is the silent file, the rest the diagnostic one)
#
# orig.tcl is made from the original's AppMain.tcl with `sed -n 'FIRST,LASTp'` for these procedures (the line numbers are
# in MECHANICS.md; keep them in this order, one after the other, then the last one):
#   gen_translation_map translate revcom triplet_transl degen_transl validation_code read_enzymes coordinate_translation
#   silent_sites do_diagnostic_site enz_site_usage_filter enz_count findenzymes, and `::tcl::mathfunc::seq_len_mod`.
# The original's own code is not part of this repository: cut it from your copy of the program.
set origfile [lindex $argv 0]
set seed [lindex $argv 1]
set nsilent [lindex $argv 2]
set ndiag [lindex $argv 3]
set ACC "/Applications/ApE.app/Contents/Resources/Scripts/Accessory Files"

# ---- stand-ins for the windowing code
proc mc {args} { if {[llength $args] > 1} { return [format [lindex $args 0] {*}[lrange $args 1 end]] }; return [lindex $args 0] }
proc tk_messageBox {args} { return ok }
proc warning_dialog {args} {}
proc add_compatibles {} {}
proc initialize_methylation {w} {}
proc toggle_methylation {w} {}
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
  set rest [lrange $args 1 end]
  for {set k 0} {$k < [llength $rest]} {incr k 2} {
    set txt [lindex $rest $k]
    set tags [lindex $rest [expr {$k+1}]]
    if {"highlight" in $tags} { append ::out "<$txt>" } else { append ::out $txt }
  }
}

source $origfile

set w .w1
set info(Accdir) $ACC
set info(clock_format) "DATE"
set info(last_analysis_xy) {0 0}
set info(translation_code) "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
set info(translation_map_lists) [gen_translation_map $info(translation_code)]
set info(Dam_Dcm_methylated) 0
set info($w,Dam_Dcm_methylated) 0
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
proc pick {n} {
  set all $::enzymes
  set chosen [list]
  while {[llength $chosen] < $n} {
    set e [lindex $all [rnd [llength $all]]]
    if {$e ni $chosen} { lappend chosen $e }
  }
  return $chosen
}
proc make_sequence {len chosen} {
  set s [randdna $len]
  foreach k {1 2 3} {
    if {[rnd 3] == 0 && [llength $chosen] > 0} {
      set e [lindex $chosen [rnd [llength $chosen]]]
      set site [concrete $::enzinfo(flatpat,$e)]
      if {[rnd 2]} { set site [revcom $site] }
      set at [rnd [expr {$len - [string length $site]}]]
      set s [string replace $s $at [expr {$at + [string length $site] - 1}] $site]
    }
  }
  if {[rnd 4] == 0} {
    set a [rnd [expr {$len - 30}]]
    set s [string replace $s $a [expr {$a + 20}] [string tolower [string range $s $a [expr {$a + 20}]]]]
  }
  if {[rnd 7] == 0} {
    set a [rnd $len]
    set s [string replace $s $a $a N]
  }
  return $s
}

proc set_case {circ sequence selection} {
  set ::seq $sequence
  set ::sel $selection
  set ::info($::w,circular) [expr {$circ ? "circular" : "linear"}]
  set ::info($::w,scanned) 0
}

proc print_case {kind k circ chosen} {
  puts "CASE $kind $k"
  puts "circular $circ"
  puts "seq $::seq"
  puts "sel [expr {$::sel eq {} ? "-" : $::sel}]"
  puts "enz [expr {[llength $chosen] == 0 ? "-" : $chosen}]"
}

proc print_sites {} {
  foreach e $::enzymes {
    if {[llength $::info($::w,$e)] > 0} { puts "sites $e [join $::info($::w,$e) { }]" }
  }
}

proc print_out {} {
  puts "out"
  foreach line [split $::out "\n"] { puts "|$line" }
  puts "END"
}

for {set k 1} {$k <= $nsilent} {incr k} {
  if {$k == 1} { set chosen {}; set len 60 } else { set chosen [pick [expr {5 + [rnd 8]}]]; set len [expr {90 + [rnd 170]}] }
  set circ [rnd 2]
  set sequence [make_sequence $len $chosen]
  set selection {}
  if {$k > 2 && [rnd 5] < 2} {
    set a [rnd [expr {$len - 30}]]
    set b [expr {$a + 12 + [rnd [expr {$len - $a - 12}]]}]
    set selection [list $a $b]
  }
  set_case $circ $sequence $selection
  findenzymes $w
  print_case silent $k $circ $chosen
  print_sites
  set ::out ""
  silent_sites $w $chosen
  print_out
}

for {set k 1} {$k <= $ndiag} {incr k} {
  set chosen [pick [expr {2 + [rnd 4]}]]
  set mode [expr {$k % 4 == 0 ? "absent" : "selected"}]
  set mm [expr {1 + [rnd 3]}]
  set len [expr {70 + [rnd 120]}]
  if {$mode eq "absent"} { set mm 1; set len [expr {50 + [rnd 50]}] }
  set circ [rnd 2]
  set sequence [make_sequence $len $chosen]
  set selection {}
  if {[rnd 5] < 2} {
    set a [rnd [expr {$len - 30}]]
    set b [expr {$a + 12 + [rnd [expr {$len - $a - 12}]]}]
    set selection [list $a $b]
  }
  set_case $circ $sequence $selection
  set info(enz_make_diag_enz) $mode
  set info(enz_make_diag_mm) $mm
  set info(enz_currently_selected) $chosen
  findenzymes $w
  puts "CASE diagnostic $k"
  puts "circular $circ"
  puts "seq $::seq"
  puts "sel [expr {$::sel eq {} ? "-" : $::sel}]"
  puts "enz [expr {[llength $chosen] == 0 ? "-" : $chosen}]"
  puts "mode $mode"
  puts "mismatches $mm"
  print_sites
  set ::out ""
  do_diagnostic_site $w
  print_out
}
