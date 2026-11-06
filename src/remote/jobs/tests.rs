use super::*;
use crate::pipeline::tools::work_folder;

fn config(dir: &Path) -> Config {
    let mut c = Config::new("t".repeat(32), dir.to_path_buf());
    c.gpus = vec![
        Gpu {
            id: 0,
            name: "GPU A".into(),
        },
        Gpu {
            id: 1,
            name: "GPU B".into(),
        },
    ];
    c.cpu_jobs = 1;
    c
}

fn store(label: &str) -> (Arc<Shared>, PathBuf) {
    let dir = work_folder(label).expect("a folder");
    (Shared::open(config(&dir), false).expect("a store"), dir)
}

#[test]
fn gpus_are_handed_out_one_job_at_a_time_lowest_first_and_come_back() {
    let mut s = Scheduler::new(&[1, 0, 1], 2);
    assert_eq!(s.take(true), Some(Slot::Gpu(0)));
    assert_eq!(s.take(true), Some(Slot::Gpu(1)));
    assert_eq!(s.take(true), None, "two GPUs, two jobs");
    assert!(!s.can_take(true));
    // CPU jobs have their own places
    assert_eq!(s.take(false), Some(Slot::Cpu));
    assert_eq!(s.take(false), Some(Slot::Cpu));
    assert_eq!(s.take(false), None);
    s.give_back(Slot::Gpu(1));
    assert_eq!(s.take(true), Some(Slot::Gpu(1)), "the one given back");
    s.give_back(Slot::Cpu);
    assert!(s.can_take(false));
    // no GPU at all: a GPU job can never start, a CPU job always has at least one place
    let none = Scheduler::new(&[], 0);
    assert!(!none.can_take(true) && none.can_take(false));
}

#[test]
fn only_ids_the_server_makes_are_valid() {
    assert!(valid_id("j6720a1b2-1f"));
    assert!(valid_id("j1-1"));
    for bad in [
        "",
        "j",
        "j-1",
        "j1-",
        "1-1",
        "J1-1",
        "j1_1",
        "j1-1/..",
        "../j1-1",
        "j1-1\n",
        "jg-1",
        "j1-1-2",
        &format!("j{}-1", "a".repeat(17)),
    ] {
        assert!(!valid_id(bad), "{bad:?}");
    }
}

#[test]
fn a_job_is_one_line_and_reads_back_whatever_its_label() {
    let job = JobInfo {
        id: "j1-2".into(),
        kind: Kind::PhyloCnn,
        label: "my\ttrees\nof \u{e9}".into(),
        status: Status::Running,
        created: 100,
        started: 105,
        finished: 0,
        gpu: Some(3),
        needs_gpu: true,
        error: "bad\tthing\r\nhappened".into(),
    };
    let line = job.to_line();
    assert!(!line.contains('\n') && line.matches('\t').count() == 9);
    let back = JobInfo::from_line(&line).expect("a job");
    assert_eq!(back.label, "my trees of \u{e9}");
    assert_eq!(back.error, "bad thing  happened");
    assert_eq!(
        (back.kind, back.status, back.gpu, back.needs_gpu),
        (Kind::PhyloCnn, Status::Running, Some(3), true)
    );
    assert_eq!(JobInfo::from_line(&back.to_line()), Some(back));
    for bad in [
        "",
        "j1-2\talign",
        "j1-2\tnope\tqueued\t1\t0\t0\t\t0\ta\tb",
        "j1-2\talign\tqueued\tx\t0\t0\t\t0\ta\tb",
    ] {
        assert!(JobInfo::from_line(bad).is_none(), "{bad:?}");
    }
}

#[test]
fn jobs_sent_in_the_same_second_are_listed_and_queued_in_the_order_they_came() {
    // the counter in an id is hexadecimal and compared as a number: as text "-10" would come before "-9"
    assert!(made_order(5, "j1-9") < made_order(5, "j1-10"));
    assert!("j1-10" < "j1-9");
    assert!(
        made_order(5, "j1-ff") < made_order(6, "j1-1"),
        "time comes first"
    );

    let (shared, dir) = store("jobs-order");
    // more than sixteen, so that the counter passes from one hex digit to two, sent within a second or two of each other
    let ids: Vec<String> = (0..20)
        .map(|n| {
            shared
                .submit(Kind::Align, &format!("job {n}"), b"x", false)
                .expect("queued")
        })
        .collect();
    let listed: Vec<String> = shared.list().into_iter().map(|j| j.id).collect();
    assert_eq!(listed, ids, "listed in the order of sending");
    // and after a restart the jobs that were waiting wait in that order too
    drop(shared);
    let again = Shared::open(config(&dir), false).expect("a store");
    let (first, slot) = again.next_job();
    assert_eq!(first, ids[0]);
    again.finish(&first, slot, Outcome::Done);
    let (second, _) = again.next_job();
    assert_eq!(second, ids[1]);
}

#[test]
fn jobs_are_queued_listed_and_removed() {
    let (shared, _dir) = store("jobs-queue");
    let a = shared
        .submit(Kind::Align, "first", b"request a", false)
        .expect("queued");
    let b = shared
        .submit(Kind::Fold, "second\nline", b"request b", true)
        .expect("queued");
    assert_ne!(a, b);
    assert!(valid_id(&a) && valid_id(&b));
    let jobs = shared.list();
    assert_eq!(
        jobs.iter().map(|j| j.id.as_str()).collect::<Vec<_>>(),
        [a.as_str(), b.as_str()]
    );
    assert!(jobs.iter().all(|j| j.status == Status::Queued));
    assert_eq!(jobs[1].label, "second line");
    assert_eq!(
        std::fs::read(shared.job_dir(&a).join("request.bin")).expect("stored"),
        b"request a"
    );
    assert_eq!(shared.overview().queued, 2);
    // a result that does not exist yet
    assert_eq!(shared.result(&a), Err(Lookup::NotReady(Status::Queued)));
    assert_eq!(shared.result("j9-9"), Err(Lookup::NotFound));
    // removing a queued job takes it out of the queue and off the disk
    shared.delete(&a).expect("removed");
    assert!(shared.get(&a).is_none() && !shared.job_dir(&a).exists());
    assert_eq!(shared.overview().queued, 1);
    assert_eq!(shared.delete(&a), Err(Lookup::NotFound));
    assert_eq!(shared.delete("../../etc"), Err(Lookup::NotFound));
}

#[test]
fn the_queue_has_a_limit() {
    let dir = work_folder("jobs-limit").expect("a folder");
    let mut c = config(&dir);
    c.max_queued = 2;
    let shared = Shared::open(c, false).expect("a store");
    shared
        .submit(Kind::Align, "1", b"x", false)
        .expect("queued");
    shared
        .submit(Kind::Align, "2", b"x", false)
        .expect("queued");
    let why = shared.submit(Kind::Align, "3", b"x", false).unwrap_err();
    assert!(why.contains("queue is full"), "{why}");
    assert_eq!(shared.list().len(), 2);
    // the refused job left nothing on the disk
    let on_disk = std::fs::read_dir(dir.join("jobs"))
        .expect("a folder")
        .count();
    assert_eq!(on_disk, 2);
}

#[test]
fn what_was_waiting_still_waits_after_a_restart_and_what_was_running_is_failed() {
    let dir = work_folder("jobs-restart").expect("a folder");
    let first = Shared::open(config(&dir), false).expect("a store");
    let a = first
        .submit(Kind::Align, "a", b"request a", false)
        .expect("queued");
    let b = first
        .submit(Kind::Tree, "b", b"request b", false)
        .expect("queued");
    let c = first
        .submit(Kind::Align, "c", b"request c", false)
        .expect("queued");
    let d = first
        .submit(Kind::Align, "d", b"request d", false)
        .expect("queued");
    // b was running, c was done with a result, d lost its request
    let mut job = first.get(&b).expect("a job");
    job.status = Status::Running;
    first.save(&job);
    let mut job = first.get(&c).expect("a job");
    job.status = Status::Done;
    job.finished = now();
    first.save(&job);
    std::fs::write(first.job_dir(&c).join("result.bin"), b"the result").expect("a result");
    std::fs::remove_file(first.job_dir(&d).join("request.bin")).expect("removed");
    drop(first);

    let second = Shared::open(config(&dir), false).expect("a store");
    let status = |id: &str| second.get(id).map(|j| j.status);
    assert_eq!(status(&a), Some(Status::Queued));
    assert_eq!(status(&b), Some(Status::Failed));
    assert_eq!(status(&c), Some(Status::Done));
    assert_eq!(status(&d), Some(Status::Failed));
    assert!(second.get(&b).expect("a job").error.contains("stopped"));
    assert_eq!(second.result(&c), Ok(b"the result".to_vec()));
    assert_eq!(second.overview().queued, 1);
    // new ids go on from where the old ones stopped
    let e = second
        .submit(Kind::Align, "e", b"x", false)
        .expect("queued");
    let count =
        |id: &str| u64::from_str_radix(id.split_once('-').expect("an id").1, 16).expect("hex");
    assert!(count(&e) > count(&d), "{e} after {d}");
    // folders that are not jobs are left alone and not read as jobs
    std::fs::create_dir_all(dir.join("jobs").join("not-a-job")).expect("a folder");
    assert_eq!(
        Shared::open(config(&dir), false)
            .expect("a store")
            .list()
            .len(),
        5
    );
}

#[test]
fn finished_jobs_are_swept_when_old_and_the_rest_are_kept() {
    let (shared, _dir) = store("jobs-sweep");
    let old = shared
        .submit(Kind::Align, "old", b"x", false)
        .expect("queued");
    let fresh = shared
        .submit(Kind::Align, "fresh", b"x", false)
        .expect("queued");
    let waiting = shared
        .submit(Kind::Align, "waiting", b"x", false)
        .expect("queued");
    let finish = |id: &str, age: u64| {
        let mut job = shared.get(id).expect("a job");
        job.status = Status::Done;
        job.finished = now() - age;
        shared.save(&job);
        shared.lock().jobs.insert(id.to_string(), job);
        shared.lock().queue.retain(|q| q != id);
    };
    let keep = shared.config.keep.as_secs();
    finish(&old, keep + 100);
    finish(&fresh, keep / 2);
    shared.sweep();
    assert!(shared.get(&old).is_none() && !shared.job_dir(&old).exists());
    assert!(shared.get(&fresh).is_some());
    assert!(
        shared.get(&waiting).is_some(),
        "a job that waits never expires"
    );
}

#[test]
fn the_next_job_goes_to_a_free_gpu_and_a_job_behind_it_may_use_the_cpu() {
    let dir = work_folder("jobs-next").expect("a folder");
    let mut c = config(&dir);
    c.gpus = vec![Gpu {
        id: 2,
        name: "only".into(),
    }];
    let shared = Shared::open(c, false).expect("a store");
    let g1 = shared
        .submit(Kind::PhyloCnn, "g1", b"x", true)
        .expect("queued");
    let g2 = shared
        .submit(Kind::PhyloCnn, "g2", b"x", true)
        .expect("queued");
    let cpu = shared
        .submit(Kind::Align, "cpu", b"x", false)
        .expect("queued");
    let (first, slot1) = shared.next_job();
    assert_eq!((first.as_str(), slot1), (g1.as_str(), Slot::Gpu(2)));
    assert_eq!(shared.get(&g1).expect("a job").gpu, Some(2));
    assert_eq!(shared.get(&g1).expect("a job").status, Status::Running);
    // the second GPU job has to wait for the GPU, so the CPU job behind it goes first
    let (second, slot2) = shared.next_job();
    assert_eq!((second.as_str(), slot2), (cpu.as_str(), Slot::Cpu));
    assert_eq!(shared.get(&g2).expect("a job").status, Status::Queued);
    let busy = shared.overview();
    assert_eq!((busy.running, busy.queued), (2, 1));
    assert_eq!(busy.gpus[0].1.as_deref(), Some(g1.as_str()));
    // when the first GPU job is over the GPU goes to the next
    shared.finish(&g1, slot1, Outcome::Failed("it broke\nbadly".into()));
    let over = shared.get(&g1).expect("a job");
    assert_eq!(
        (over.status, over.error.as_str()),
        (Status::Failed, "it broke badly")
    );
    let (third, slot3) = shared.next_job();
    assert_eq!((third.as_str(), slot3), (g2.as_str(), Slot::Gpu(2)));
}

#[test]
fn stopping_a_running_job_leaves_it_cancelled_and_removing_it_again_forgets_it() {
    let (shared, _dir) = store("jobs-stop");
    let id = shared
        .submit(Kind::Align, "x", b"x", false)
        .expect("queued");
    let (_, slot) = shared.next_job();
    shared.delete(&id).expect("asked to stop");
    // nothing is removed while the process runs; the worker sees the request and ends the job
    assert_eq!(shared.get(&id).expect("a job").status, Status::Running);
    shared.finish(&id, slot, Outcome::Failed("killed".into()));
    assert_eq!(shared.get(&id).expect("a job").status, Status::Cancelled);
    assert!(shared.job_dir(&id).exists());
    shared.delete(&id).expect("removed");
    assert!(shared.get(&id).is_none() && !shared.job_dir(&id).exists());
    // a job that finished its work just as it was stopped keeps its result
    let id = shared
        .submit(Kind::Align, "y", b"x", false)
        .expect("queued");
    let (_, slot) = shared.next_job();
    shared.delete(&id).expect("asked to stop");
    shared.finish(&id, slot, Outcome::Done);
    assert_eq!(shared.get(&id).expect("a job").status, Status::Done);
}

#[test]
fn what_was_sent_is_deleted_when_the_job_is_over_and_only_so_many_jobs_are_kept() {
    let dir = work_folder("jobs-kept").expect("a folder");
    let mut c = config(&dir);
    c.max_kept = 3;
    let shared = Shared::open(c, false).expect("a store");
    let a = shared
        .submit(Kind::Align, "a", b"request a", false)
        .expect("queued");
    assert!(shared.job_dir(&a).join("request.bin").is_file());
    let (_, slot) = shared.next_job();
    shared.finish(&a, slot, Outcome::Done);
    assert!(
        !shared.job_dir(&a).join("request.bin").exists(),
        "a finished job does not keep what was sent"
    );
    // three jobs are kept (one done, two waiting): the fourth is refused before anything is written
    shared
        .submit(Kind::Align, "b", b"x", false)
        .expect("queued");
    shared
        .submit(Kind::Align, "c", b"x", false)
        .expect("queued");
    let why = shared
        .submit(Kind::Align, "d", b"x", false)
        .expect_err("no room");
    assert!(why.contains("keeps 3 jobs"), "{why}");
    let folders = std::fs::read_dir(dir.join("jobs")).expect("jobs").count();
    assert_eq!(folders, 3, "the refused job left nothing behind");
    // removing a finished job makes room
    shared.delete(&a).expect("removed");
    shared
        .submit(Kind::Align, "e", b"x", false)
        .expect("queued");
}

// a folder can only be opened (to set its age) like this on unix
#[cfg(unix)]
#[test]
fn a_half_made_job_folder_is_cleared_away_when_old_and_its_number_is_never_reused() {
    let dir = work_folder("jobs-orphan").expect("a folder");
    let jobs = dir.join("jobs");
    let (old, fresh) = (jobs.join("j1-ff"), jobs.join("j1-100"));
    for folder in [&old, &fresh] {
        std::fs::create_dir_all(folder).expect("a folder");
    }
    // eight days old: older than the seven that a result is kept
    let long_ago = SystemTime::now() - Duration::from_secs(8 * 24 * 3600);
    File::open(&old)
        .and_then(|f| f.set_modified(long_ago))
        .expect("an old folder");
    let shared = Shared::open(config(&dir), false).expect("a store");
    assert!(!old.exists(), "the old one is cleared away");
    assert!(fresh.exists(), "a new one may be a job being made");
    // a new job takes a number above every folder's, readable or not
    let id = shared
        .submit(Kind::Align, "next", b"x", false)
        .expect("queued");
    assert!(counter_of(&id) > 0x100, "{id}");
}

#[test]
fn a_job_that_finished_just_before_the_server_stopped_keeps_its_result() {
    let dir = work_folder("jobs-just-done").expect("a folder");
    let shared = Shared::open(config(&dir), false).expect("a store");
    let done = shared
        .submit(Kind::Align, "done", b"x", false)
        .expect("queued");
    let lost = shared
        .submit(Kind::Align, "lost", b"x", false)
        .expect("queued");
    // the first job is over, and the second is running
    let (first, slot) = shared.next_job();
    assert_eq!(first, done);
    shared.finish(&first, slot, Outcome::Done);
    let (running, _) = shared.next_job();
    assert_eq!(running, lost);
    // the second job's process wrote a whole result (by a rename) and the server stopped before it could say so
    std::fs::write(shared.job_dir(&lost).join("result.bin"), b"a result").expect("a result");
    drop(shared);
    let again = Shared::open(config(&dir), false).expect("a store");
    assert_eq!(again.get(&done).expect("a job").status, Status::Done);
    let job = again.get(&lost).expect("a job");
    assert_eq!(job.status, Status::Done, "its result was there");
    assert_eq!(again.result(&lost).expect("a result"), b"a result");
}

#[cfg(unix)]
#[test]
fn a_folder_that_was_already_there_keeps_its_permissions_and_new_ones_are_closed_to_others() {
    use std::os::unix::fs::PermissionsExt;
    let mode = |p: &Path| std::fs::metadata(p).expect("a folder").permissions().mode() & 0o777;
    let dir = work_folder("jobs-modes").expect("a folder");
    std::fs::set_permissions(&dir, std::fs::Permissions::from_mode(0o755)).expect("permissions");
    let shared = Shared::open(config(&dir), false).expect("a store");
    // the folder the person chose is left as it was (it may be shared); what the server made is for the owner only
    assert_eq!(mode(&dir), 0o755);
    assert_eq!(mode(&dir.join("jobs")), 0o700);
    let id = shared
        .submit(Kind::Align, "x", b"x", false)
        .expect("queued");
    assert_eq!(mode(&shared.job_dir(&id)), 0o700);
}
