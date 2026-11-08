# Runs the original's Text Map code (`do_text_map` and everything it calls), unchanged, against REAL Tk text widgets, and prints
# the cases and what the original wrote: its text and the ranges of every tag. orig.tcl is cut out of AppMain.tcl with
# cut_procs.tcl. This is how src/seq/textmap/golden.txt was made (MECHANICS.md, section 8, "Text Map").
#
# The system Tcl (8.5.9) loads its own Tk, so the sequence window's text widget (rows of equal width, feature tags with the
# `<<Metadata>>` and `<<Revcolors>>` bindings, as `put_tags` makes them) and the result window's text widget are the real thing;
# only the windowing code around them (analysis window, geometry, menus) is stood in for.
#
# usage: tclsh cut_procs.tcl AppMain.tcl revcom validation_code read_enzymes not_in_list initialize_methylation \
#          toggle_methylation findenzymes ix2bp bp2ix end_ix textarea_get get_tags_exact get_cds linear_enz_sort translate \
#          gen_translation_map uppercase_cdna copy_tags do_text_map ::tcl::mathfunc::seq_len_mod > orig.tcl
#        tclsh text_map_driver.tcl orig.tcl <seed> <cases> > golden.txt
set origfile [lindex $argv 0]
set seed [lindex $argv 1]
set ncases [lindex $argv 2]
set ACC "/Applications/ApE.app/Contents/Resources/Scripts/Accessory Files"

package require Tk
wm withdraw .
fconfigure stdout -buffering line

# ---- stand-ins for the windowing code (the text widgets themselves are real)
proc mc {args} { if {[llength $args] > 1} { return [format [lindex $args 0] {*}[lrange $args 1 end]] }; return [lindex $args 0] }
proc sputs {args} {}
proc tk_messageBox {args} { return ok }
proc warning_dialog {args} {}
proc add_compatibles {} {}
rename wm ::tkwm
proc wm {args} {
  switch -- [lindex $args 0] { title - withdraw { return [tkwm {*}$args] } }
  return ""
}
proc window_on_screen {args} {}
proc focus {args} {}
proc textbox_geometry {args} { return "" }
proc add_output_menu {args} { return "" }
proc update_windows_menu {} {}
set ::analysis 0
proc new_analysis_window {w title name} {
  incr ::analysis
  toplevel .an$::analysis
  tkwm withdraw .an$::analysis
  return .an$::analysis
}
proc create_textframe {a width} {
  frame $a.textframe
  text $a.textframe.text -wrap none -width 40 -height 10
  pack $a.textframe $a.textframe.text
}

source $origfile

set w .w1
set info(Accdir) $ACC
set info(clock_format) "DATE"
set info(last_analysis_xy) {0 0}
set info(translation_code) FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
set info(translation_map_lists) [gen_translation_map $info(translation_code)]
set dialogblock 0
set ok 0
read_enzymes [file join $ACC Enzymes Default_Enzymes.txt]

expr {srand($seed)}
proc rnd {n} { return [expr {int(rand()*$n)}] }
proc pick {list} { return [lindex $list [rnd [llength $list]]] }
proc randdna {n} { set s ""; for {set i 0} {$i < $n} {incr i} { append s [string index ACGT [rnd 4]] }; return $s }
proc concrete {pat} {
  set out ""
  foreach c [split $pat ""] {
    set bases [string map {A A C C G G T T R AG Y CT S CG W AT K GT M AC B CGT D AGT H ACT V ACG N ACGT} $c]
    append out [string index $bases [rnd [string length $bases]]]
  }
  return $out
}
# a random sequence with the sites of several enzymes put in (on either strand, some twice)
proc make_sequence {len} {
  set s [randdna $len]
  set n [expr {2 + [rnd 10]}]
  for {set k 0} {$k < $n} {incr k} {
    set e [lindex $::enzymes [rnd [llength $::enzymes]]]
    set site [concrete $::enzinfo(flatpat,$e)]
    if {[rnd 2]} { set site [revcom $site] }
    if {[string length $site] >= $len} { continue }
    set at [rnd [expr {$len - [string length $site]}]]
    set s [string replace $s $at [expr {$at + [string length $site] - 1}] $site]
  }
  return $s
}
# DNA with exons in capitals and introns (gt...ag) in small letters, a small-letter lead and a small-letter tail after a stop,
# so that "Uppercase as Genes" has something to do
proc make_gene_like {} {
  set s [string tolower [randdna [rnd 12]]]
  set exons [expr {1 + [rnd 3]}]
  for {set k 0} {$k < $exons} {incr k} {
    append s [randdna [expr {6 + [rnd 40]}]]
    if {$k < $exons - 1} { append s "gt[string tolower [randdna [rnd 15]]]ag" }
  }
  append s [pick {TGA TAA TAG}]
  append s [string tolower [randdna [rnd 12]]]
  return $s
}
# some ambiguity codes and some small letters
proc sprinkle {s} {
  set n [string length $s]
  for {set k 0} {$k < 1 + [rnd 4]} {incr k} {
    set at [rnd $n]
    set s [string replace $s $at $at [pick {N R Y S W K M B D H V}]]
  }
  if {[rnd 2]} {
    set at [rnd $n]
    set to [expr {$at + [rnd 30]}]
    set s [string replace $s $at $to [string tolower [string range $s $at $to]]]
  }
  return $s
}

# ---- the sequence window: a real text widget in rows of equal width, tags as `put_tags` makes them
proc open_document {seq title width circ features sel} {
  catch {destroy .w1}
  toplevel .w1
  tkwm withdraw .w1
  tkwm title .w1 $title
  text .w1.textarea
  set rows [list]
  for {set i 0} {$i < [string length $seq]} {incr i $width} { lappend rows [string range $seq $i [expr {$i + $width - 1}]] }
  .w1.textarea insert 1.0 [join $rows "\n"]
  set len [string length $seq]
  set t .w1.textarea
  foreach f $features {
    lassign $f i dir hidden fwd rev type ranges fn name
    set tag "f$i#"
    foreach {s e} $ranges { add_range $t $tag $s $e }
    $t tag configure $tag -background [expr {$hidden ? "" : ($dir ? $rev : $fwd)}] -foreground {} -font {} -underline {}
    $t tag bind $tag <<Revcolors>> [list $fwd $rev]
    $t tag bind $tag <<Metadata>> [list $name {} $type $dir {} {}]
  }
  # `put_tags` (as a file is opened, L2029) makes the tags of the joined features after all the others
  foreach f $features {
    lassign $f i dir hidden fwd rev type ranges fn name
    if {$fn ne ""} {
      set fntag "fn$i#"
      foreach {s e} $fn { add_range $t $fntag $s $e }
      $t tag configure $fntag -background {} -foreground {} -font {} -underline {}
    }
  }
  if {$sel ne ""} {
    foreach {s e} $sel { add_range $t sel $s $e }
  }
}
proc add_range {t tag s e} {
  if {[$t compare [bp2ix $t $e] > [end_ix $t]]} {
    $t tag add $tag [bp2ix $t $s] [end_ix $t]
  } else {
    $t tag add $tag [bp2ix $t $s] [bp2ix $t $e]
  }
}

# ---- random features: ranges as the GenBank reader makes them (start-1, end), and the `fn` range of a joined one
proc random_segments {from to count} {
  # `count` increasing, non-empty segments inside [from, to); sometimes two touch
  set room [expr {$to - $from}]
  if {$room < 2 * $count} { set count 1 }
  if {$room < 2} { return {} }
  set cuts [list]
  while {[llength $cuts] < 2 * $count} {
    set c [expr {$from + [rnd [expr {$room + 1}]]}]
    lappend cuts $c
  }
  set cuts [lsort -integer $cuts]
  set segs [list]
  foreach {a b} $cuts {
    if {$a == $b} { continue }
    lappend segs $a $b
  }
  if {[llength $segs] == 0} { return [list $from $to] }
  return $segs
}
proc random_feature {i len circ} {
  set dir [rnd 2]
  set hidden [expr {[rnd 5] == 0}]
  set fwd [format "#%02x%02x%02x" [rnd 256] [rnd 256] [rnd 256]]
  set rev [format "#%02x%02x%02x" [rnd 256] [rnd 256] [rnd 256]]
  set type [pick {CDS CDS CDS gene promoter misc_feature rep_origin primer_bind terminator cds}]
  set name [pick {amp lacZ "T7 promoter" ori "my gene" a "long feature name that goes on and on and on" "x-y" "b c d" "(odd) name" "CDS 1"}]
  set fn ""
  set kind [rnd 8]
  if {$circ && $len > 12 && $kind == 0} {
    # one range through the origin, as one location: {a len 0 b}
    set a [expr {$len / 2 + [rnd [expr {$len - $len / 2 - 1}]]}]
    set b [expr {1 + [rnd [expr {$len / 2}]]}]
    set ranges [list $a $len 0 $b]
  } elseif {$circ && $len > 24 && $kind == 1} {
    # joined, through the origin: the parts before it, then the parts after it
    set pre [random_segments [expr {$len / 2}] $len [expr {1 + [rnd 2]}]]
    set post [random_segments 0 [expr {$len / 2}] [expr {1 + [rnd 2]}]]
    if {[llength $pre] == 0 || [llength $post] == 0} { return [random_feature $i $len 0] }
    set ranges [concat $pre $post]
    set fn [fn_of $ranges $len $circ]
  } else {
    set count [expr {[rnd 3] == 0 ? 1 + [rnd 3] : 1}]
    set ranges [random_segments 0 $len $count]
    set fn [fn_of $ranges $len $circ]
  }
  return [list $i $dir $hidden $fwd $rev $type $ranges $fn $name]
}
# the range of the `fn` tag of a feature of several parts, as `feature_table_to_taglist` makes it (L41228-41245): from the first
# base to the last, or, when the first part begins after the last one ends on a circle, from the first base to the end and
# from the start to the last; two parts that are one location through the origin have none
proc fn_of {ranges len circ} {
  if {[llength $ranges] <= 2} { return "" }
  if {$circ && [lindex $ranges 0] > [lindex $ranges end]} {
    if {[llength $ranges] == 4 && [lindex $ranges 1] == $len && [lindex $ranges 2] == 0} { return "" }
    return [list [lindex $ranges 0] $len 0 [lindex $ranges end]]
  }
  return [list [lindex $ranges 0] [lindex $ranges end]]
}

proc rand_char {} { return [pick [list " " + * : | . ^ -]] }

set frames {index enzymes dna translation sec_strand genes features space}
proc shuffled {list} {
  set out [list]
  while {[llength $list]} {
    set k [rnd [llength $list]]
    lappend out [lindex $list $k]
    set list [lreplace $list $k $k]
  }
  return $out
}

proc set_info {key value} { set ::info($key) $value }

for {set k 1} {$k <= $ncases} {incr k} {
  # progress on stderr; a run that says nothing on stderr has been seen to lose all its output on stdout (a Tk quirk)
  if {$k % 100 == 1} { puts stderr "text_map_driver: case $k of $ncases" }
  set len [expr {1 + [rnd 260]}]
  if {[rnd 12] == 0} { set len [expr {1 + [rnd 6]}] }
  set circ [rnd 2]
  set methylated [rnd 2]
  if {[rnd 3] == 0 && $len > 30} {
    set seq [make_gene_like]
    set len [string length $seq]
  } else {
    set seq [make_sequence $len]
    if {[rnd 4] == 0} { set seq [sprinkle $seq] }
  }
  set width [expr {5 + [rnd 70]}]
  set title [pick {"pUC19.gb" "a" "a long name for a sequence document that is much longer than one line of the text map"}]
  set features [list]
  set nfeat [expr {[rnd 3] == 0 ? 0 : [rnd 8]}]
  for {set i 0} {$i < $nfeat} {incr i} { lappend features [random_feature $i $len $circ] }
  # the part looked at
  set region [expr {[rnd 3] == 0 ? "all" : "selection"}]
  set sel ""
  if {$region eq "selection" || [rnd 4] == 0} {
    set a [rnd $len]
    set b [expr {$a + 1 + [rnd [expr {$len - $a}]]}]
    set sel [list $a $b]
  }
  set picked [expr {[rnd 2] ? $::enzymes : [lrange [shuffled $::enzymes] 0 [rnd 12]]}]
  # the settings the dialog holds
  set info(enz_text_width) [pick {10 12 20 33 50 60 100 100 100 150}]
  set info(enz_text_space) [expr {1 + [rnd 7]}]
  set info(enz_text_show_space) [rnd 2]
  set info(enz_text_dna_numbers) [pick {None Left Right Both}]
  set info(enz_text_copy_highlight) [rnd 2]
  set info(enz_text_2nd) [rnd 2]
  set info(enz_text_show_trans) [rnd 2]
  set info(enz_text_trans) [pick {1_1 3_1 6_1 1_3 3_3 6_3}]
  # the original's site finder starts on the first five bases and has not much to look at in a sequence shorter than a dozen,
  # which is the business of its own tests (List Enzymes), not of this one
  set info(enz_text_show_enz) [expr {$len >= 12 && [rnd 3] > 0}]
  set info(enz_text_genes) [rnd 2]
  set info(enz_text_graphic_map) 0
  set info(enz_text_show_index) [rnd 2]
  set info(enz_text_index_char1) [rand_char]
  set info(enz_text_index_char2) [rand_char]
  set info(enz_text_index_spacing1) [expr {1 + [rnd 25]}]
  set info(enz_text_index_spacing2) [expr {1 + [rnd 12]}]
  set info(enz_text_index_nums) [rnd 2]
  set info(enz_text_show_features) [expr {[rnd 4] > 0}]
  set info(enz_text_show_hidden_features) [rnd 2]
  set info(enz_text_crop_features) [rnd 2]
  set info(enz_text_dialog_order) [shuffled $frames]
  # the second strand is made from the DNA line before it: put before the DNA it stops with an error (as in case 1 of the
  # kind), which is kept in only one case in ten
  if {[rnd 10] > 0} {
    set o $info(enz_text_dialog_order)
    set a [lsearch $o sec_strand]
    set b [lsearch $o dna]
    if {$a < $b} { lset o $a dna; lset o $b sec_strand; set info(enz_text_dialog_order) $o }
  }
  set info($w,circular) [expr {$circ ? "circular" : "linear"}]
  set info($w,scanned) 0
  set info($w,Dam_Dcm_methylated) $methylated
  open_document $seq $title $width $circ $features $sel
  findenzymes $w
  puts "CASE textmap $k"
  puts "circular $circ"
  puts "seq $seq"
  puts "title $title"
  puts "region [expr {$region eq "all" || $sel eq "" ? "all" : $sel}]"
  puts "rowwidth $width"
  puts "selection [expr {$sel eq "" ? "-" : $sel}]"
  puts "methylated $methylated"
  puts "enz [join $picked { }]"
  puts "settings width=$info(enz_text_width) space=$info(enz_text_space) show_space=$info(enz_text_show_space) numbers=$info(enz_text_dna_numbers) copy=$info(enz_text_copy_highlight) second=$info(enz_text_2nd) show_trans=$info(enz_text_show_trans) trans=$info(enz_text_trans) show_enz=$info(enz_text_show_enz) genes=$info(enz_text_genes) show_index=$info(enz_text_show_index) char1=[string map {{ } _} $info(enz_text_index_char1)] char2=[string map {{ } _} $info(enz_text_index_char2)] sp1=$info(enz_text_index_spacing1) sp2=$info(enz_text_index_spacing2) nums=$info(enz_text_index_nums) show_features=$info(enz_text_show_features) hidden=$info(enz_text_show_hidden_features) crop=$info(enz_text_crop_features)"
  puts "order [join $info(enz_text_dialog_order) { }]"
  foreach f $features {
    lassign $f i dir hidden fwd rev type ranges fn name
    puts "feature $i $dir $hidden $fwd $rev $type | [join $ranges { }] | [join $fn { }] | $name"
  }
  foreach e $::enzymes {
    if {[llength $info($w,$e)] > 0} { puts "sites $e [join $info($w,$e) { }]" }
  }
  set analyze [expr {$sel eq "" ? "all" : "selection"}]
  if {$region eq "all"} { set analyze all }
  if {[catch {set tw [do_text_map $w $picked $analyze]} err]} {
    puts "ERROR $err"
  } else {
    puts "out"
    foreach line [split [$tw get 1.0 end-1c] "\n"] { puts "|$line" }
    puts "tags"
    foreach tag [lsort [$tw tag names]] {
      if {$tag eq "highlite" || $tag eq "sel"} { continue }
      set r [$tw tag ranges $tag]
      if {[llength $r] > 0} {
        set bg [$tw tag cget $tag -background]
        puts "tag $tag [expr {$bg eq "" ? "-" : $bg}] [join $r { }]"
      }
    }
  }
  puts "END"
  catch {destroy .an$::analysis}
}
flush stdout
exit
