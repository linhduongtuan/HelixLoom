//! Hybrid assembly after Unicycler (Wick, Judd, Gorrie and Holt, *PLoS Comput Biol* 2017, 13:e1005595): an assembly graph from short reads, as
//! GFA, is joined up with long reads. Written from the paper alone (Unicycler is GPL licensed; none of its source was read or copied, and its
//! README could not be read): where the paper gives no rule the program has its own, and `docs/HYBRID.md` says which. What is here:
//!
//! * `input`: the graph as GFA, with the depth of each segment;
//! * `copies`: how many copies of the genome's sequence each segment is: single copy, or repeat (from depth and from the links);
//! * `route`: the paths through the repeats of the graph that join two single-copy segments;
//! * `bridges`: where long reads join two single-copy segments, found with the program's mapper;
//! * `quality`: a score for a bridge (the paper lists what it depends on, not the formula; this one is the program's own);
//! * `apply`: the bridges, best first, put into the graph, and the graph made compact again;
//! * `rotate`: a circular replicon turned so that it starts at a start gene;
//! * `assemble`, `output`, `cli`: the steps in order, the files, and `helixloom hybrid`.

mod apply;
pub mod assemble;
mod bridges;
pub mod cli;
mod copies;
pub mod input;
pub mod output;
mod quality;
pub mod rotate;
mod route;
pub use apply::Mode;

#[cfg(test)]
pub(crate) mod sim;
#[cfg(test)]
mod tests;
