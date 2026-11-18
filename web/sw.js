// The service worker of HelixLoom's page (docs/WEB.md, "Installed, and without a network").
//
// It keeps the files of the build it came with, so that the page opens when there is no network: the program runs in the
// browser, and everything but the work that is sent to a workstation goes on. scripts/build-web.sh fills in the two
// lines below: the version is a checksum of the build, and the files are the page and what the page needs.
//
// What it answers is only what is listed here. Every other request is left to the browser as if there were no worker:
// above all the calls to a workstation, which carry a token and bring back the person's results, and are never kept.
const VERSION = "__VERSION__";
const FILES = [__FILES__];

const SCOPE = self.registration.scope;
// The caches of this page, whatever their build: a site may have another HelixLoom page in another folder, whose
// caches are that page's and are named after its folder.
const OF_THIS_PAGE = "helixloom " + SCOPE + " ";
const KEPT = OF_THIS_PAGE + VERSION;
const OURS = new Set(FILES.map((file) => new URL(file, SCOPE).href));
// the page, under the two names it is asked for by
const PAGES = new Set([SCOPE, new URL("index.html", SCOPE).href]);

/** An answer as one of its own: a browser refuses to show a page from an answer that came by way of a redirect. */
function whole(answer) {
  if (!answer.redirected) return answer;
  return new Response(answer.body, { status: answer.status, statusText: answer.statusText, headers: answer.headers });
}

/** Keep every file of the build, or none of it: a page must never be put together from the files of two builds. */
async function keepTheBuild() {
  // asked of the server itself, not of what the browser has kept by its own rules (a server that says nothing about
  // keeping would otherwise hand the page of the build before to a worker of this one)
  const answers = await Promise.all(FILES.map((file) => fetch(new Request(file, { cache: "reload" }))));
  const missing = answers.find((answer) => !answer.ok);
  if (missing) throw new Error("not kept: " + missing.url + " (" + missing.status + ")");
  const cache = await caches.open(KEPT);
  await Promise.all(FILES.map((file, at) => cache.put(new URL(file, SCOPE).href, whole(answers[at]))));
}

self.addEventListener("install", (event) => {
  event.waitUntil(
    keepTheBuild().catch((why) =>
      // what was kept of a build that could not be kept whole goes again
      caches.delete(KEPT).then(() => {
        throw why;
      }),
    ),
  );
});

// The caches of the builds before this one go when this one takes over; a cache that is not this page's is not touched.
self.addEventListener("activate", (event) => {
  event.waitUntil(
    caches.keys().then((names) => {
      const old = names.filter((name) => name.startsWith(OF_THIS_PAGE) && name !== KEPT);
      return Promise.all(old.map((name) => caches.delete(name)));
    }),
  );
});

/** The kept copy of `address`, from this build's cache and no other. */
function kept(address) {
  return caches.open(KEPT).then((cache) => cache.match(address));
}

self.addEventListener("fetch", (event) => {
  const request = event.request;
  if (request.method !== "GET") return;
  if (request.mode === "navigate") {
    // the page is asked for with what it is to show (?view=trees): that is no part of which page it is
    const page = new URL(request.url);
    page.search = "";
    page.hash = "";
    if (!PAGES.has(page.href)) return;
    // the newest page when the server answers (a new build is seen at once), the kept one when it does not; and when
    // there is neither, what a browser shows when it cannot reach a page, as if there were no worker
    event.respondWith(
      fetch(request).catch(() => kept(SCOPE).then((page) => page || Response.error())),
    );
    return;
  }
  if (!OURS.has(request.url)) return;
  // a file of this build is what it was when it was kept: its name is a checksum of what is in it
  event.respondWith(kept(request.url).then((found) => found || fetch(request)));
});
