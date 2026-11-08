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
    assert_eq!(s.take(1), Some(Slot::Gpus(vec![0])));
    assert_eq!(s.take(1), Some(Slot::Gpus(vec![1])));
    assert_eq!(s.take(1), None, "two GPUs, two jobs");
    assert!(!s.can_take(1));
    // CPU jobs have their own places
    assert_eq!(s.take(0), Some(Slot::Cpu));
    assert_eq!(s.take(0), Some(Slot::Cpu));
    assert_eq!(s.take(0), None);
    s.give_back(&Slot::Gpus(vec![1]));
    assert_eq!(s.take(1), Some(Slot::Gpus(vec![1])), "the one given back");
    s.give_back(&Slot::Cpu);
    assert!(s.can_take(0));
    // no GPU at all: a GPU job can never start, a CPU job always has at least one place
    let none = Scheduler::new(&[], 0);
    assert!(!none.can_take(1) && none.can_take(0));
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
        gpus: vec![3],
        gpus_wanted: 1,
        error: "bad\tthing\r\nhappened".into(),
    };
    let line = job.to_line();
    assert!(!line.contains('\n') && line.matches('\t').count() == 9);
    let back = JobInfo::from_line(&line).expect("a job");
    assert_eq!(back.label, "my trees of \u{e9}");
    assert_eq!(back.error, "bad thing  happened");
    assert_eq!(
        (back.kind, back.status, back.gpus.clone(), back.gpus_wanted),
        (Kind::PhyloCnn, Status::Running, vec![3], 1)
    );
    assert!(back.needs_gpu());
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
fn a_job_with_several_gpus_is_one_line_too_and_a_line_from_before_that_still_reads() {
    let mut job = JobInfo {
        id: "j1-3".into(),
        kind: Kind::PhyloCnn,
        label: "many trees".into(),
        status: Status::Running,
        created: 1,
        started: 2,
        finished: 0,
        gpus: vec![0, 1, 3],
        gpus_wanted: 3,
        error: String::new(),
    };
    let line = job.to_line();
    assert_eq!(line.split('\t').nth(6), Some("0,1,3"));
    assert_eq!(line.split('\t').nth(7), Some("3"));
    assert_eq!(JobInfo::from_line(&line), Some(job.clone()));
    // a job that waits for two GPUs has none yet and wants two
    job.status = Status::Queued;
    job.gpus.clear();
    job.gpus_wanted = 2;
    assert_eq!(JobInfo::from_line(&job.to_line()), Some(job));
    // a line as a server wrote it when a job had one GPU at most: a number, and a 1 or a 0
    let old =
        JobInfo::from_line("j1-4\tphylocnn\trunning\t1\t2\t0\t2\t1\tx\t").expect("an old line");
    assert_eq!((old.gpus, old.gpus_wanted), (vec![2], 1));
    let none = JobInfo::from_line("j1-5\talign\tqueued\t1\t0\t0\t\t0\tx\t").expect("an old line");
    assert_eq!((none.gpus, none.gpus_wanted), (vec![], 0));
    // numbers nobody could mean are not a job
    for bad in [
        "j1-6\talign\tqueued\t1\t0\t0\t64\t1\tx\t",
        "j1-6\talign\tqueued\t1\t0\t0\t0,\t1\tx\t",
        "j1-6\talign\tqueued\t1\t0\t0\ta\t1\tx\t",
        "j1-6\talign\tqueued\t1\t0\t0\t\t65\tx\t",
        "j1-6\talign\tqueued\t1\t0\t0\t\t-1\tx\t",
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
                .submit(Kind::Align, &format!("job {n}"), b"x", 0)
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
    again.finish(&first, &slot, Outcome::Done);
    let (second, _) = again.next_job();
    assert_eq!(second, ids[1]);
}

#[test]
fn jobs_are_queued_listed_and_removed() {
    let (shared, _dir) = store("jobs-queue");
    let a = shared
        .submit(Kind::Align, "first", b"request a", 0)
        .expect("queued");
    let b = shared
        .submit(Kind::Fold, "second\nline", b"request b", 1)
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
    shared.submit(Kind::Align, "1", b"x", 0).expect("queued");
    shared.submit(Kind::Align, "2", b"x", 0).expect("queued");
    let why = shared.submit(Kind::Align, "3", b"x", 0).unwrap_err();
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
        .submit(Kind::Align, "a", b"request a", 0)
        .expect("queued");
    let b = first
        .submit(Kind::Tree, "b", b"request b", 0)
        .expect("queued");
    let c = first
        .submit(Kind::Align, "c", b"request c", 0)
        .expect("queued");
    let d = first
        .submit(Kind::Align, "d", b"request d", 0)
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
    let e = second.submit(Kind::Align, "e", b"x", 0).expect("queued");
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
    let old = shared.submit(Kind::Align, "old", b"x", 0).expect("queued");
    let fresh = shared
        .submit(Kind::Align, "fresh", b"x", 0)
        .expect("queued");
    let waiting = shared
        .submit(Kind::Align, "waiting", b"x", 0)
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
        .submit(Kind::PhyloCnn, "g1", b"x", 1)
        .expect("queued");
    let g2 = shared
        .submit(Kind::PhyloCnn, "g2", b"x", 1)
        .expect("queued");
    let cpu = shared.submit(Kind::Align, "cpu", b"x", 0).expect("queued");
    let (first, slot1) = shared.next_job();
    assert_eq!(
        (first.as_str(), slot1.clone()),
        (g1.as_str(), Slot::Gpus(vec![2]))
    );
    assert_eq!(shared.get(&g1).expect("a job").gpus, vec![2]);
    assert_eq!(shared.get(&g1).expect("a job").status, Status::Running);
    // the second GPU job has to wait for the GPU, so the CPU job behind it goes first
    let (second, slot2) = shared.next_job();
    assert_eq!((second.as_str(), slot2), (cpu.as_str(), Slot::Cpu));
    assert_eq!(shared.get(&g2).expect("a job").status, Status::Queued);
    let busy = shared.overview();
    assert_eq!((busy.running, busy.queued), (2, 1));
    assert_eq!(busy.gpus[0].1.as_deref(), Some(g1.as_str()));
    // when the first GPU job is over the GPU goes to the next
    shared.finish(&g1, &slot1, Outcome::Failed("it broke\nbadly".into()));
    let over = shared.get(&g1).expect("a job");
    assert_eq!(
        (over.status, over.error.as_str()),
        (Status::Failed, "it broke badly")
    );
    let (third, slot3) = shared.next_job();
    assert_eq!((third.as_str(), slot3), (g2.as_str(), Slot::Gpus(vec![2])));
}

/// A store of a workstation with four GPUs, numbered 0 to 3, and room for one job without a GPU.
fn four_gpus(label: &str) -> Arc<Shared> {
    let dir = work_folder(label).expect("a folder");
    let mut c = config(&dir);
    c.gpus = (0..4)
        .map(|id| Gpu {
            id,
            name: format!("GPU {id}"),
        })
        .collect();
    Shared::open(c, false).expect("a store")
}

#[test]
fn a_job_that_wants_several_gpus_is_given_them_all_at_once_and_the_status_page_shows_each() {
    let shared = four_gpus("jobs-several");
    let big = shared
        .submit(Kind::PhyloCnn, "big", b"x", 3)
        .expect("queued");
    let small = shared
        .submit(Kind::PhyloCnn, "small", b"x", 1)
        .expect("queued");
    let (first, slot) = shared.next_job();
    assert_eq!(
        (first.as_str(), slot.clone()),
        (big.as_str(), Slot::Gpus(vec![0, 1, 2]))
    );
    assert_eq!(shared.get(&big).expect("a job").gpus, vec![0, 1, 2]);
    // the job's process is shown exactly these
    assert_eq!(slot.visible_devices(), "0,1,2");
    let (second, slot2) = shared.next_job();
    assert_eq!(
        (second.as_str(), slot2.clone()),
        (small.as_str(), Slot::Gpus(vec![3]))
    );
    let busy = shared.overview();
    let on = |gpu: usize| busy.gpus[gpu].1.clone();
    assert_eq!(on(0).as_deref(), Some(big.as_str()));
    assert_eq!(on(1).as_deref(), Some(big.as_str()));
    assert_eq!(on(2).as_deref(), Some(big.as_str()));
    assert_eq!(on(3).as_deref(), Some(small.as_str()));
    // all of them come back, and are used again
    shared.finish(&big, &slot, Outcome::Done);
    let again = shared
        .submit(Kind::PhyloCnn, "again", b"x", 2)
        .expect("queued");
    let (third, slot3) = shared.next_job();
    assert_eq!(
        (third.as_str(), slot3),
        (again.as_str(), Slot::Gpus(vec![0, 1]))
    );
    shared.finish(&small, &slot2, Outcome::Done);
}

#[test]
fn a_job_that_waits_for_several_gpus_is_not_passed_by_jobs_that_want_one() {
    let shared = four_gpus("jobs-order-gpus");
    // GPU 0 and GPU 1 are in use
    let busy1 = shared
        .submit(Kind::PhyloCnn, "b1", b"x", 2)
        .expect("queued");
    let (_, busy_slot) = shared.next_job();
    // a job wants three and only two are free; a job that wants one comes after it
    let wide = shared
        .submit(Kind::PhyloCnn, "wide", b"x", 3)
        .expect("queued");
    let narrow = shared
        .submit(Kind::PhyloCnn, "narrow", b"x", 1)
        .expect("queued");
    // a job with no GPU at all is not held up by either
    let cpu = shared.submit(Kind::Align, "cpu", b"x", 0).expect("queued");
    assert!(JobInfo::needs_gpu(&shared.get(&wide).expect("a job")));
    let (next, slot) = shared.next_job();
    assert_eq!(
        (next.as_str(), slot.clone()),
        (cpu.as_str(), Slot::Cpu),
        "the CPU job goes"
    );
    assert_eq!(
        shared.get(&narrow).expect("a job").status,
        Status::Queued,
        "the narrow job waits for the wide one"
    );
    assert_eq!(shared.get(&wide).expect("a job").status, Status::Queued);
    // when the first job is over, the wide one starts, and the narrow one after it, on what is left
    shared.finish(&busy1, &busy_slot, Outcome::Done);
    let (next, slot) = shared.next_job();
    assert_eq!(
        (next.as_str(), slot),
        (wide.as_str(), Slot::Gpus(vec![0, 1, 2]))
    );
    let (next, slot) = shared.next_job();
    assert_eq!(
        (next.as_str(), slot),
        (narrow.as_str(), Slot::Gpus(vec![3]))
    );
}

#[test]
fn removing_the_job_at_the_head_of_the_line_lets_the_jobs_behind_it_start_at_once() {
    let shared = four_gpus("jobs-order-removed");
    // GPU 0 and GPU 1 are in use
    shared
        .submit(Kind::PhyloCnn, "busy", b"x", 2)
        .expect("queued");
    let (_, _busy_slot) = shared.next_job();
    // the job at the head wants three GPUs and two are free; a job that wants one is behind it
    let wide = shared
        .submit(Kind::PhyloCnn, "wide", b"x", 3)
        .expect("queued");
    let narrow = shared
        .submit(Kind::PhyloCnn, "narrow", b"x", 1)
        .expect("queued");
    // a worker waits for a job to start, as the server's workers do
    let (send, receive) = std::sync::mpsc::channel();
    let worker = Arc::clone(&shared);
    std::thread::spawn(move || {
        let _ = send.send(worker.next_job());
    });
    assert!(
        receive.recv_timeout(Duration::from_millis(300)).is_err(),
        "nothing starts while the wide job is first"
    );
    // the person removes the wide job: the worker must not go on sleeping while a GPU is free
    shared.delete(&wide).expect("removed");
    let (next, slot) = receive
        .recv_timeout(Duration::from_secs(10))
        .expect("the narrow job starts as soon as the wide one is gone");
    assert_eq!(
        (next.as_str(), slot),
        (narrow.as_str(), Slot::Gpus(vec![2]))
    );
}

#[test]
fn a_job_queued_for_more_gpus_than_the_server_now_has_gets_what_there_is() {
    let shared = four_gpus("jobs-fewer");
    let id = shared.submit(Kind::PhyloCnn, "x", b"x", 8).expect("queued");
    let (first, slot) = shared.next_job();
    assert_eq!(
        (first.as_str(), slot),
        (id.as_str(), Slot::Gpus(vec![0, 1, 2, 3]))
    );
    assert_eq!(shared.get(&id).expect("a job").gpus_wanted, 8);
}

#[test]
fn stopping_a_running_job_leaves_it_cancelled_and_removing_it_again_forgets_it() {
    let (shared, _dir) = store("jobs-stop");
    let id = shared.submit(Kind::Align, "x", b"x", 0).expect("queued");
    let (_, slot) = shared.next_job();
    shared.delete(&id).expect("asked to stop");
    // nothing is removed while the process runs; the worker sees the request and ends the job
    assert_eq!(shared.get(&id).expect("a job").status, Status::Running);
    shared.finish(&id, &slot, Outcome::Failed("killed".into()));
    assert_eq!(shared.get(&id).expect("a job").status, Status::Cancelled);
    assert!(shared.job_dir(&id).exists());
    shared.delete(&id).expect("removed");
    assert!(shared.get(&id).is_none() && !shared.job_dir(&id).exists());
    // a job that finished its work just as it was stopped keeps its result
    let id = shared.submit(Kind::Align, "y", b"x", 0).expect("queued");
    let (_, slot) = shared.next_job();
    shared.delete(&id).expect("asked to stop");
    shared.finish(&id, &slot, Outcome::Done);
    assert_eq!(shared.get(&id).expect("a job").status, Status::Done);
}

#[test]
fn what_was_sent_is_deleted_when_the_job_is_over_and_only_so_many_jobs_are_kept() {
    let dir = work_folder("jobs-kept").expect("a folder");
    let mut c = config(&dir);
    c.max_kept = 3;
    let shared = Shared::open(c, false).expect("a store");
    let a = shared
        .submit(Kind::Align, "a", b"request a", 0)
        .expect("queued");
    assert!(shared.job_dir(&a).join("request.bin").is_file());
    let (_, slot) = shared.next_job();
    shared.finish(&a, &slot, Outcome::Done);
    assert!(
        !shared.job_dir(&a).join("request.bin").exists(),
        "a finished job does not keep what was sent"
    );
    // three jobs are kept (one done, two waiting): the fourth is refused before anything is written
    shared.submit(Kind::Align, "b", b"x", 0).expect("queued");
    shared.submit(Kind::Align, "c", b"x", 0).expect("queued");
    let why = shared
        .submit(Kind::Align, "d", b"x", 0)
        .expect_err("no room");
    assert!(why.contains("keeps 3 jobs"), "{why}");
    let folders = std::fs::read_dir(dir.join("jobs")).expect("jobs").count();
    assert_eq!(folders, 3, "the refused job left nothing behind");
    // removing a finished job makes room
    shared.delete(&a).expect("removed");
    shared.submit(Kind::Align, "e", b"x", 0).expect("queued");
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
    let id = shared.submit(Kind::Align, "next", b"x", 0).expect("queued");
    assert!(counter_of(&id) > 0x100, "{id}");
}

#[test]
fn a_job_that_finished_just_before_the_server_stopped_keeps_its_result() {
    let dir = work_folder("jobs-just-done").expect("a folder");
    let shared = Shared::open(config(&dir), false).expect("a store");
    let done = shared.submit(Kind::Align, "done", b"x", 0).expect("queued");
    let lost = shared.submit(Kind::Align, "lost", b"x", 0).expect("queued");
    // the first job is over, and the second is running
    let (first, slot) = shared.next_job();
    assert_eq!(first, done);
    shared.finish(&first, &slot, Outcome::Done);
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
    let id = shared.submit(Kind::Align, "x", b"x", 0).expect("queued");
    assert_eq!(mode(&shared.job_dir(&id)), 0o700);
}
