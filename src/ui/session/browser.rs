//! The browser's storage for the page's site (`localStorage`) as the table a store keeps things in (`medium.rs`): in a page
//! there is no settings folder, and without this the web version would keep nothing between visits. What is kept is under
//! names of its own, so that another page of the same site keeps its things apart. The browser gives a site a few megabytes:
//! the appearance, the places of the windows and the open sequences of a plasmid editor fit many times over, and a set of
//! sequences that does not fit is told, not half kept (`session.rs` writes the list of what is open last).

/// What is said when the table cannot take a text.
#[cfg_attr(not(any(target_arch = "wasm32", test)), allow(dead_code))]
pub const FULL: &str =
    "The browser's storage for this site is full or may not be written, so your work is not being \
                        kept: save what you need (File menu).";

/// Every name of this program in the browser's table begins with this.
const OWN: &str = "helixloom/";

/// The name of the table's entry that keeps `name`.
#[cfg_attr(not(any(target_arch = "wasm32", test)), allow(dead_code))]
pub fn entry(name: &str) -> String {
    format!("{OWN}{name}")
}

/// The name that an entry of the table keeps, if it is one of this program's.
#[cfg_attr(not(any(target_arch = "wasm32", test)), allow(dead_code))]
pub fn kept_name(entry: &str) -> Option<&str> {
    entry.strip_prefix(OWN).filter(|name| !name.is_empty())
}

#[cfg(target_arch = "wasm32")]
mod local {
    use super::super::medium::Shelf;
    use super::{entry, kept_name, FULL};

    /// One of the two tables a browser gives a page. It is looked up each time: a page whose storage is switched off has
    /// none, and says so on a write.
    #[derive(Debug, Clone, Copy)]
    pub enum Web {
        /// The site's (`localStorage`): kept until the person clears it. The session and the preferences are here.
        Site,
        /// The tab's (`sessionStorage`): the tab's alone, there when the page is loaded again in it, and gone when the
        /// tab is closed. The workstation's address and token are here, and nowhere else.
        Tab,
    }

    impl Web {
        fn table(self) -> Option<web_sys::Storage> {
            let window = web_sys::window()?;
            match self {
                Web::Site => window.local_storage(),
                Web::Tab => window.session_storage(),
            }
            .ok()?
        }

        /// This table as a shelf, if the page has it.
        fn shelf(self) -> Option<std::rc::Rc<dyn Shelf>> {
            self.table()
                .map(|_| std::rc::Rc::new(self) as std::rc::Rc<dyn Shelf>)
        }
    }

    impl Shelf for Web {
        fn get(&self, name: &str) -> Option<String> {
            self.table()?.get_item(&entry(name)).ok()?
        }

        fn set(&self, name: &str, text: &str) -> Result<(), String> {
            self.table()
                .ok_or_else(|| FULL.to_string())?
                .set_item(&entry(name), text)
                .map_err(|_| FULL.to_string())
        }

        fn remove(&self, name: &str) {
            if let Some(table) = self.table() {
                let _ = table.remove_item(&entry(name));
            }
        }

        fn names(&self) -> Vec<String> {
            let Some(table) = self.table() else {
                return Vec::new();
            };
            let count = table.length().unwrap_or(0);
            (0..count)
                .filter_map(|i| table.key(i).ok().flatten())
                .filter_map(|key| kept_name(&key).map(str::to_string))
                .collect()
        }
    }

    /// The table of the page's site, if the page has one.
    pub fn shelf() -> Option<std::rc::Rc<dyn Shelf>> {
        Web::Site.shelf()
    }

    /// The table of this tab, if the page has one.
    pub fn tab_shelf() -> Option<std::rc::Rc<dyn Shelf>> {
        Web::Tab.shelf()
    }
}

/// The browser's database (IndexedDB) as the larger table of `large.rs`: the results of the pages are kept in it, and nothing else.
/// It is opened, and read whole, before the page starts (`open`), because a store reads at once and a database answers later; a
/// page whose browser has none, or cannot open it, goes on without, and keeps its results in the small table as it did.
#[cfg(target_arch = "wasm32")]
mod database {
    use super::super::large::{Large, Op, Persist, ROOM};
    use js_sys::{Array, Function, Promise};
    use std::cell::RefCell;
    use std::rc::Rc;
    use wasm_bindgen::prelude::*;
    use wasm_bindgen::JsCast;
    use wasm_bindgen_futures::JsFuture;
    use web_sys::{IdbDatabase, IdbRequest, IdbTransactionMode};

    /// The database of this program, and its one table of names and texts.
    const NAME: &str = "helixloom";
    const TABLE: &str = "results";
    const VERSION: u32 = 1;
    /// How long the opening, and the reading of what is in it, are waited for before the page starts without it.
    const OPEN_PATIENCE_MS: i32 = 4_000;
    const READ_PATIENCE_MS: i32 = 20_000;

    type Done = Box<dyn FnOnce(Result<(), String>)>;
    type Finish = Rc<dyn Fn(Result<(), String>)>;

    thread_local! {
        /// The table that `open` made, for the stores that the page makes.
        static LARGE: RefCell<Option<Large>> = const { RefCell::new(None) };
    }

    /// What went wrong, in words: the name and the message of a DOM exception, or what a script threw.
    fn describe(error: &JsValue) -> String {
        match error.dyn_ref::<web_sys::DomException>() {
            Some(e) => format!("{}: {}", e.name(), e.message()),
            None => error
                .as_string()
                .unwrap_or_else(|| "the database stopped".to_string()),
        }
    }

    /// A promise that is kept with what a request gave, and broken with why it failed.
    fn promise_of(request: &IdbRequest) -> Promise {
        let request = request.clone();
        Promise::new(&mut move |resolve: Function, reject: Function| {
            let read = request.clone();
            let ok = Closure::once_into_js(move |_: web_sys::Event| {
                let result = read.result().unwrap_or(JsValue::UNDEFINED);
                let _ = resolve.call1(&JsValue::NULL, &result);
            });
            let failed = request.clone();
            let bad = Closure::once_into_js(move |_: web_sys::Event| {
                let why = match failed.error() {
                    Ok(Some(e)) => format!("{}: {}", e.name(), e.message()),
                    _ => "the database stopped".to_string(),
                };
                let _ = reject.call1(&JsValue::NULL, &JsValue::from_str(&why));
            });
            request.set_onsuccess(Some(ok.unchecked_ref()));
            request.set_onerror(Some(bad.unchecked_ref()));
        })
    }

    /// `promise`, or a failure when it takes more than `ms` milliseconds: a database that does not answer must not keep the page
    /// from starting.
    async fn within(promise: Promise, ms: i32) -> Result<JsValue, String> {
        let window = web_sys::window().ok_or("there is no browser window")?;
        let too_long = Promise::new(&mut |_resolve: Function, reject: Function| {
            let late = Closure::once_into_js(move || {
                let why = JsValue::from_str("the database did not answer in time");
                let _ = reject.call1(&JsValue::NULL, &why);
            });
            let _ = window
                .set_timeout_with_callback_and_timeout_and_arguments_0(late.unchecked_ref(), ms);
        });
        let race = Promise::race(&Array::of2(&promise, &too_long));
        JsFuture::from(race).await.map_err(|why| describe(&why))
    }

    /// The database, opened (and made, the first time).
    async fn opened() -> Result<IdbDatabase, String> {
        let window = web_sys::window().ok_or("there is no browser window")?;
        let factory = window
            .indexed_db()
            .map_err(|why| describe(&why))?
            .ok_or("this browser has no database for pages")?;
        let request = factory
            .open_with_u32(NAME, VERSION)
            .map_err(|why| describe(&why))?;
        let making = request.clone();
        let upgrade = Closure::once_into_js(move |_: web_sys::Event| {
            if let Ok(db) = making.result().and_then(|r| r.dyn_into::<IdbDatabase>()) {
                let _ = db.create_object_store(TABLE);
            }
        });
        request.set_onupgradeneeded(Some(upgrade.unchecked_ref()));
        let db = within(promise_of(&request), OPEN_PATIENCE_MS).await?;
        db.dyn_into::<IdbDatabase>()
            .map_err(|_| "the database did not open".to_string())
    }

    /// Everything the table holds, each name with its text.
    async fn read(db: &IdbDatabase) -> Result<Vec<(String, String)>, String> {
        let table = db
            .transaction_with_str_and_mode(TABLE, IdbTransactionMode::Readonly)
            .and_then(|transaction| transaction.object_store(TABLE))
            .map_err(|why| describe(&why))?;
        // both are asked in one transaction before anything is waited for, so that they are of the same moment
        let keys = table.get_all_keys().map_err(|why| describe(&why))?;
        let values = table.get_all().map_err(|why| describe(&why))?;
        let both = Promise::all(&Array::of2(&promise_of(&keys), &promise_of(&values)));
        let both = Array::from(&within(both, READ_PATIENCE_MS).await?);
        let (keys, values) = (Array::from(&both.get(0)), Array::from(&both.get(1)));
        Ok(keys
            .iter()
            .zip(values.iter())
            .filter_map(|(key, value)| Some((key.as_string()?, value.as_string()?)))
            .collect())
    }

    /// The database as what keeps the changes of a table (`large.rs`).
    #[derive(Debug)]
    struct Idb {
        db: IdbDatabase,
    }

    impl Idb {
        /// Ask for one change in a transaction of its own, and have `finish` called once with how it went (a change that is
        /// refused takes no other with it).
        fn start(&self, op: &Op, finish: &Finish) -> Result<(), String> {
            let transaction = self
                .db
                .transaction_with_str_and_mode(TABLE, IdbTransactionMode::Readwrite)
                .map_err(|why| describe(&why))?;
            let asked = (|| -> Result<(), JsValue> {
                let table = transaction.object_store(TABLE)?;
                match op {
                    Op::Put(name, text) => {
                        table.put_with_key(&JsValue::from_str(text), &JsValue::from_str(name))?
                    }
                    Op::Delete(name) => table.delete(&JsValue::from_str(name))?,
                    Op::Clear => table.clear()?,
                };
                Ok(())
            })();
            if let Err(why) = asked {
                let _ = transaction.abort();
                return Err(describe(&why));
            }
            let (complete, stopped) = (finish.clone(), finish.clone());
            let ok = Closure::once_into_js(move |_: web_sys::Event| complete(Ok(())));
            let seen = transaction.clone();
            // an abort is how a full database says so: the transaction's error names it
            let bad = Closure::once_into_js(move |_: web_sys::Event| {
                let why = seen.error().map_or_else(
                    || "the database stopped".to_string(),
                    |e| format!("{}: {}", e.name(), e.message()),
                );
                stopped(Err(why));
            });
            transaction.set_oncomplete(Some(ok.unchecked_ref()));
            transaction.set_onabort(Some(bad.unchecked_ref()));
            // a closure that is not called (the other of the two) is not freed: a few bytes for each save
            Ok(())
        }
    }

    impl Persist for Idb {
        fn apply(&self, op: Op, done: Done) {
            let done = Rc::new(RefCell::new(Some(done)));
            let finish: Finish = Rc::new(move |outcome| {
                if let Some(done) = done.borrow_mut().take() {
                    done(outcome);
                }
            });
            if let Err(why) = self.start(&op, &finish) {
                finish(Err(why));
            }
        }
    }

    /// Open the database and read what it holds into the larger table, for the stores the page makes. Without a database, or
    /// when it does not open in time, the page goes on without it.
    pub async fn open() {
        // a page whose browser keeps no `localStorage` has no store at all (`shelf`): the database would be read for nothing
        if super::local::shelf().is_none() {
            return;
        }
        let table = async {
            let db = opened().await?;
            let loaded = read(&db).await?;
            Ok::<Large, String>(Large::new(Rc::new(Idb { db }), loaded, ROOM))
        };
        match table.await {
            Ok(table) => LARGE.with(|large| *large.borrow_mut() = Some(table)),
            Err(why) => log::warn!(
                "HelixLoom cannot use the browser's database, so the results of the pages are kept in the small table: {why}"
            ),
        }
    }

    /// The larger table, if the page has one.
    pub fn installed() -> Option<Large> {
        LARGE.with(|large| large.borrow().clone())
    }

    /// Have the page drawn again when the database refuses a change: nothing else would draw a frame, and the person is told on one.
    pub fn wake_on_refusal(ctx: &eframe::egui::Context) {
        if let Some(large) = installed() {
            let ctx = ctx.clone();
            large.set_waker(Rc::new(move || ctx.request_repaint()));
        }
    }
}

/// The table for the page's site: the small one, with the larger one for the results when the page has it.
#[cfg(target_arch = "wasm32")]
pub fn shelf() -> Option<std::rc::Rc<dyn super::medium::Shelf>> {
    let small = local::shelf()?;
    Some(match database::installed() {
        Some(large) => std::rc::Rc::new(super::large::Roomy::new(small, Some(large))),
        None => small,
    })
}

#[cfg(target_arch = "wasm32")]
pub use database::{open as open_database, wake_on_refusal};
#[cfg(target_arch = "wasm32")]
pub use local::tab_shelf;
