//! Who may use what: the GPUs, handed out to the jobs that need them (one or several to a job, the lowest numbers first), and
//! a limited number of jobs that use none at once.

/// What a worker may use for a job.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Slot {
    Cpu,
    /// These GPUs, by their numbers, lowest first; there is at least one.
    Gpus(Vec<usize>),
}

impl Slot {
    /// The GPUs a job is shown (`CUDA_VISIBLE_DEVICES`): their numbers with commas between, nothing for a job with none.
    pub fn visible_devices(&self) -> String {
        match self {
            Slot::Cpu => String::new(),
            Slot::Gpus(ids) => ids
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(","),
        }
    }

    /// The numbers of the GPUs in this slot (none for a CPU place).
    pub fn gpus(&self) -> &[usize] {
        match self {
            Slot::Cpu => &[],
            Slot::Gpus(ids) => ids,
        }
    }
}

/// Hands out the GPUs to the jobs, and a limited number of CPU jobs at once.
#[derive(Debug)]
pub struct Scheduler {
    free_gpus: Vec<usize>,
    /// How many GPUs the server has in all: no job can have more.
    total_gpus: usize,
    free_cpu: usize,
}

impl Scheduler {
    pub fn new(gpus: &[usize], cpu_jobs: usize) -> Scheduler {
        let mut free_gpus = gpus.to_vec();
        free_gpus.sort_unstable();
        free_gpus.dedup();
        Scheduler {
            total_gpus: free_gpus.len(),
            free_gpus,
            free_cpu: cpu_jobs.max(1),
        }
    }

    /// How many GPUs a job that wants `wanted` of them is given: what it wants, but no more than there are (a job that was
    /// queued when the server had more GPUs must not wait for ones that are gone). `wanted` 0 is a job with no GPU.
    fn needed(&self, wanted: usize) -> usize {
        if wanted == 0 {
            0
        } else {
            wanted.min(self.total_gpus.max(1))
        }
    }

    /// Whether a job that wants `wanted` GPUs (0: none) could start now.
    pub fn can_take(&self, wanted: usize) -> bool {
        match self.needed(wanted) {
            0 => self.free_cpu > 0,
            n => self.free_gpus.len() >= n,
        }
    }

    /// The places for a job that wants `wanted` GPUs (0: none): the lowest free GPUs, or one of the CPU places.
    pub fn take(&mut self, wanted: usize) -> Option<Slot> {
        if !self.can_take(wanted) {
            return None;
        }
        match self.needed(wanted) {
            0 => {
                self.free_cpu -= 1;
                Some(Slot::Cpu)
            }
            n => Some(Slot::Gpus(self.free_gpus.drain(..n).collect())),
        }
    }

    pub fn give_back(&mut self, slot: &Slot) {
        match slot {
            Slot::Cpu => self.free_cpu += 1,
            Slot::Gpus(ids) => {
                self.free_gpus.extend(ids);
                self.free_gpus.sort_unstable();
            }
        }
    }
}

/// How many GPUs a job gets: what its sender asked for (0 means every GPU of the server), but at least one and no more than
/// the server has or the job can use (`most`: a job of three trees has no use for eight GPUs).
pub fn gpus_for(asked: usize, have: usize, most: usize) -> usize {
    let wanted = if asked == 0 { have } else { asked.min(have) };
    wanted.min(most).max(1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_job_gets_the_lowest_free_gpus_and_they_come_back_in_order() {
        let mut s = Scheduler::new(&[3, 1, 2, 0, 1], 2);
        assert_eq!(s.take(2), Some(Slot::Gpus(vec![0, 1])));
        assert_eq!(s.take(1), Some(Slot::Gpus(vec![2])));
        assert!(!s.can_take(2), "only GPU 3 is free");
        assert_eq!(s.take(2), None);
        s.give_back(&Slot::Gpus(vec![0, 1]));
        // the ones that came back are the lowest again, and GPU 3 is still free after them
        assert_eq!(s.take(3), Some(Slot::Gpus(vec![0, 1, 3])));
        assert_eq!(s.take(1), None, "every GPU is busy now");
        s.give_back(&Slot::Gpus(vec![2]));
        assert_eq!(s.take(1), Some(Slot::Gpus(vec![2])));
    }

    #[test]
    fn jobs_without_a_gpu_take_the_cpu_places_and_never_a_gpu() {
        let mut s = Scheduler::new(&[0], 2);
        assert_eq!(s.take(0), Some(Slot::Cpu));
        assert_eq!(s.take(0), Some(Slot::Cpu));
        assert!(!s.can_take(0));
        assert!(s.can_take(1), "the GPU is free all the same");
        s.give_back(&Slot::Cpu);
        assert!(s.can_take(0));
        // at least one CPU place whatever was asked
        assert!(Scheduler::new(&[], 0).can_take(0));
    }

    #[test]
    fn a_job_never_waits_for_more_gpus_than_the_server_has() {
        let mut s = Scheduler::new(&[4, 5], 1);
        // it was queued when there were four
        assert!(s.can_take(4));
        assert_eq!(s.take(4), Some(Slot::Gpus(vec![4, 5])));
        assert!(!s.can_take(4), "none left");
        // a server with no GPU gives none, whatever is wanted
        assert!(!Scheduler::new(&[], 1).can_take(2));
    }

    #[test]
    fn a_job_is_shown_the_gpus_it_was_given_by_number() {
        assert_eq!(Slot::Cpu.visible_devices(), "");
        assert_eq!(Slot::Gpus(vec![2]).visible_devices(), "2");
        assert_eq!(Slot::Gpus(vec![0, 1, 3]).visible_devices(), "0,1,3");
        assert!(Slot::Cpu.gpus().is_empty());
        assert_eq!(Slot::Gpus(vec![1, 4]).gpus(), [1, 4]);
    }

    #[test]
    fn the_number_of_gpus_a_job_gets_is_what_it_asked_for_within_what_there_is_and_what_it_can_use()
    {
        // (asked, the server has, the job can use) -> given
        for (asked, have, most, given) in [
            (1, 4, 100, 1),
            (2, 4, 100, 2),
            (9, 4, 100, 4),
            (0, 4, 100, 4),
            (0, 4, 3, 3),
            (4, 4, 1, 1),
            (0, 0, 100, 1),
            (3, 2, 0, 1),
        ] {
            assert_eq!(
                gpus_for(asked, have, most),
                given,
                "{asked} of {have}, {most}"
            );
        }
    }
}
