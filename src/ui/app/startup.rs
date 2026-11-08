//! What `--view NAME` shows: each name opens a window or a page with its example, so that every screen can be pictured.

use super::*;

impl App {
    /// Open the window named on the command line (used to capture each screen).
    pub(super) fn apply_startup_view(&mut self, view: &str) {
        if let Some(ws) = Workspace::from_view(view) {
            self.load_example(ws);
            return;
        }
        match view {
            "enzymes" => self.enzyme_win.open = true,
            "digest" => {
                self.selected = ["EcoRI", "HindIII", "KpnI"]
                    .iter()
                    .map(|s| s.to_string())
                    .collect();
                self.digest_win.open = true;
            }
            "find" => {
                self.find_win.open = true;
                self.find_win.text = "GAATTC".into();
                self.find_highlight_all();
            }
            "enzlist" => self.run_quick_list(Quick::SiteUsage),
            "enzlist-dialog" => self.open_list_dialog(),
            "linear" => self.map_linear = true,
            "assembly" | "assembly-before" => {
                self.tab = Tab::Assembly;
                self.assembly.load_demo();
                if view == "assembly" {
                    self.assembly.run();
                }
            }
            "gibson" => {
                self.tab = Tab::Assembly;
                self.assembly.load_demo();
                self.assembly.set_preset(Preset::Gibson);
            }
            "orfs" => {
                self.prefs.orf_min = 300;
                self.prefs.orf_strand = Strand::Both;
                self.orf_win.list_open = true;
            }
            "translate" => {
                self.docs[0].select(326, 446);
                self.docs[0].trans = TransMode::Forward;
                self.open_translate();
            }
            "features" => {
                self.docs[0].select(150, 170);
                self.open_feature_window(None);
            }
            "pcr" => {
                self.docs[0].select(140, 520);
                self.pcr_win.open = true;
                self.pcr_win.fwd = "TAATACGACTCACTATAGGGATGGGAGCAG".into();
                // a reverse primer: the reverse complement of bases 481 to 506, so a product forms
                let reverse = revcom(self.docs[0].doc.slice(480, 506));
                self.pcr_win.rev = String::from_utf8_lossy(&reverse).into_owned();
                self.find_products();
            }
            "silent" | "diagnostic" => {
                // the coding stretch of the demo plasmid, and enzymes that make a few places
                self.docs[0].select(326, 446);
                self.selected = [
                    "BamHI", "EcoRV", "HindIII", "KpnI", "NheI", "PstI", "SacI", "XhoI",
                ]
                .iter()
                .map(|s| s.to_string())
                .collect();
                if view == "silent" {
                    self.run_silent_sites();
                } else {
                    self.site_win.diagnostic_open = true;
                    self.site_win.mismatches = 2;
                    self.run_diagnostic_sites();
                }
            }
            "calc" => {
                self.calc_win.open = true;
                self.calc_win.text = "ATGGCGTTTGAACGGCTCAA".into();
            }
            "calculator" => {
                self.calculator_win.open = true;
                self.calculator_win.text = "(1 + sqrt(5)) / 2".into();
            }
            "molcalc" | "molcalc-mix" | "molcalc-molarity" => {
                self.molcalc_win.open = true;
                self.molcalc_win.load_example();
                self.molcalc_win
                    .show_tab(view.strip_prefix("molcalc-").unwrap_or(""));
            }
            "tree" | "tree-ml" | "tree-bayes" => {
                self.docs = crate::sample::tree_demo()
                    .into_iter()
                    .map(DocView::new)
                    .collect();
                self.tree_win.open = true;
                match view {
                    "tree-ml" => self.tree_win.analysis = Analysis::Likelihood,
                    "tree-bayes" => self.tree_win.analysis = Analysis::Bayesian,
                    _ => {}
                }
                self.build_tree();
            }
            "phylocnn" => {
                self.phylocnn_win.open = true;
                self.load_phylocnn_example();
                self.run_phylocnn();
            }
            "phylocnn-many" => {
                self.phylocnn_win.open = true;
                self.load_phylocnn_example();
                self.phylocnn_win.newick = crate::sample::tree_sample(24);
                self.run_phylocnn();
            }
            "species" => {
                self.species_win.open = true;
                self.load_species_example();
                self.run_species();
            }
            "treetools" => {
                self.treetools_win.open = true;
                self.load_treetools_demo();
            }
            "primer3" => {
                self.p3_win.open = true;
                self.docs[0].select(224, 280); // the MCS feature of the demo plasmid
                self.run_primer3();
            }
            "reads" => {
                self.show_reads(&crate::sample::demo_fastq(), "demo reads");
                self.clean_reads();
            }
            "translation-rows" => self.docs[0].trans = TransMode::Six,
            "about" => self.about_open = true,
            "signals-reverse" => {
                self.load_example(Workspace::Signals);
                self.platform.signals.reverse = true;
            }
            "trees-circular" => {
                self.load_example(Workspace::Trees);
                let o = &mut self.platform.trees.options;
                o.layout = crate::seq::treelayout::Layout::Circular;
                o.groups = true;
            }
            "palette" => self.palette.toggle(),
            "shortcuts" => self.shortcuts_open = true,
            _ => {}
        }
    }
}
