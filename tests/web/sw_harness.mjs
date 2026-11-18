// The service worker of the web version (web/sw.js, as scripts/build-web.sh fills it in), run in a stand-in for a browser:
// a cache storage in memory, a network that can be switched off, and the three events a browser sends a worker. The
// page's own short script, which asks for the worker, is run too.
// Run by tests/web.rs:  node sw_harness.mjs DIST/sw.js DIST/manifest.webmanifest DIST/index.html
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import vm from "node:vm";

const [workerPath, manifestPath, pagePath] = process.argv.slice(2);
const source = readFileSync(workerPath, "utf8");
// a page that is not at the root of its site: every address has to be worked out from where the worker is
const SCOPE = "http://pages.test/helixloom/";

/** A request as a worker makes one: an address, and how the browser's own keeping is to be used. */
class Request {
  constructor(address, how = {}) {
    this.url = new URL(address, SCOPE).href;
    this.cache = how.cache;
  }
}

/** An answer a worker makes itself. */
class Response {
  constructor(body, how = {}) {
    Object.assign(this, { body, ok: true, redirected: false, made: true, ...how });
  }
  static error() {
    return { ok: false, type: "error", from: "nowhere" };
  }
}

/** A browser for one worker: what it keeps, what it asked the network for and how, and whether there is a network. */
function browser() {
  const world = {
    stores: new Map(),
    asked: [],
    online: true,
    missing: new Set(),
    redirected: new Set(),
    // addresses the browser has no room to keep
    noRoomFor: new Set(),
    listeners: {},
  };
  const network = async (request) => {
    const url = typeof request === "string" ? request : request.url;
    world.asked.push({ url, cache: request.cache });
    if (!world.online) throw new TypeError("no network");
    const ok = !world.missing.has(url);
    return { ok, status: ok ? 200 : 404, url, body: "the bytes of " + url, redirected: world.redirected.has(url), from: "network" };
  };
  class Cache {
    constructor() {
      this.held = new Map();
    }
    async put(address, answer) {
      assert.equal(typeof address, "string", "a kept answer is kept under an address");
      assert.ok(!answer.redirected, "an answer that came by a redirect is not kept as it came: " + address);
      if (world.noRoomFor.has(address)) throw new Error("no room to keep " + address);
      this.held.set(address, { ...answer, from: "kept" });
    }
    // as a browser's: the whole address has to be the same, what is asked after the "?" too
    async match(request) {
      const url = typeof request === "string" ? new URL(request, SCOPE).href : request.url;
      return this.held.get(url);
    }
  }
  const caches = {
    open: async (name) => {
      if (!world.stores.has(name)) world.stores.set(name, new Cache());
      return world.stores.get(name);
    },
    keys: async () => [...world.stores.keys()],
    delete: async (name) => world.stores.delete(name),
    // a worker that asks every cache at once could be answered from another build's: it must name its own
    match: async () => assert.fail("the worker asked all the caches, not its own"),
  };
  const self = {
    location: new URL("sw.js", SCOPE),
    registration: { scope: SCOPE },
    addEventListener: (type, listener) => {
      world.listeners[type] = listener;
    },
  };
  vm.runInNewContext(source, { self, caches, fetch: network, URL, Promise, Set, Error, Request, Response, console });
  /** Send the worker an event that it may ask to be waited for. */
  world.send = async (type) => {
    let waited = Promise.resolve();
    world.listeners[type]({ waitUntil: (promise) => (waited = promise) });
    await waited;
  };
  /** Ask for an address as a page would; what the worker answered with, or "left" when it left it to the browser. */
  world.ask = async (url, how = {}) => {
    let answer = "left";
    const request = { url, method: how.method ?? "GET", mode: how.mode ?? "no-cors" };
    world.listeners.fetch({ request, respondWith: (promise) => (answer = promise) });
    return answer === "left" ? "left" : await answer;
  };
  world.kept = () => [...world.stores.values()].flatMap((cache) => [...cache.held.keys()]).sort();
  return world;
}

const held = [];
async function holds(what, check) {
  await check();
  held.push(what);
}

const listed = /const FILES = \[([^\]]*)\]/.exec(source);
assert.ok(listed, "the worker lists the files of its build");
const FILES = listed[1].split(",").map((f) => f.trim().replace(/^"|"$/g, "")).filter(Boolean);
const urlOf = (file) => new URL(file, SCOPE).href;
const script = FILES.find((f) => f.endsWith(".js"));
const program = FILES.find((f) => f.endsWith(".wasm"));
assert.ok(script && program && FILES.includes("./"), "the page, its script and its program are listed: " + FILES);

await holds("it keeps every file of its build, asked of the server itself, in a cache of its page and its build", async () => {
  const b = browser();
  await b.send("install");
  assert.equal(b.stores.size, 1);
  const [name] = b.stores.keys();
  assert.match(name, /^helixloom http:\/\/pages\.test\/helixloom\/ [0-9a-f]{16}$/);
  assert.deepEqual(b.kept(), FILES.map(urlOf).sort());
  assert.deepEqual(b.asked.map((a) => a.url).sort(), FILES.map(urlOf).sort());
  assert.ok(b.asked.every((a) => a.cache === "reload"), "not what the browser kept by its own rules");
});

await holds("a build of which one file cannot be had is not kept at all", async () => {
  for (const lost of [urlOf(program), urlOf(script), SCOPE]) {
    const b = browser();
    b.missing.add(lost);
    await assert.rejects(b.send("install"), lost);
    assert.deepEqual(b.kept(), [], lost);
    assert.equal(b.stores.size, 0, "and no empty cache is left: " + lost);
  }
});

await holds("a build that the browser has no room to keep whole leaves nothing kept", async () => {
  const b = browser();
  b.noRoomFor.add(urlOf(program));
  await assert.rejects(b.send("install"), /no room/);
  assert.equal(b.stores.size, 0, "no cache with half a build in it: " + b.kept());
});

await holds("a page that the server sends by way of a redirect is kept as an answer of its own", async () => {
  const b = browser();
  b.redirected.add(SCOPE);
  await b.send("install");
  b.online = false;
  const page = await b.ask(SCOPE, { mode: "navigate" });
  assert.ok(page.made && !page.redirected && page.from === "kept");
  assert.equal(page.body, "the bytes of " + SCOPE);
});

await holds("with a network the page is the newest; without one it is the kept one", async () => {
  const b = browser();
  await b.send("install");
  for (const page of [SCOPE, SCOPE + "?view=trees", SCOPE + "index.html", SCOPE + "index.html?dark", SCOPE + "#top"]) {
    b.online = true;
    assert.equal((await b.ask(page, { mode: "navigate" })).from, "network", page);
    b.online = false;
    const kept = await b.ask(page, { mode: "navigate" });
    assert.equal(kept.from, "kept", page);
    assert.equal(kept.url, SCOPE, page);
  }
});

await holds("without a network and with nothing kept the page is what a browser shows when it reaches nothing", async () => {
  const b = browser();
  b.online = false;
  const answer = await b.ask(SCOPE, { mode: "navigate" });
  assert.ok(answer && answer.type === "error", "an answer, and one that says there is none");
});

await holds("a file of the build comes from what is kept, and the network is not asked", async () => {
  const b = browser();
  await b.send("install");
  b.asked.length = 0;
  b.online = false;
  for (const file of FILES.filter((f) => f !== "./")) {
    assert.equal((await b.ask(urlOf(file))).from, "kept", file);
  }
  assert.deepEqual(b.asked, []);
});

await holds("a file that is kept by no cache yet is fetched", async () => {
  // a worker that is running although its cache is gone (the person cleared the site's data)
  const b = browser();
  assert.equal((await b.ask(urlOf(script))).from, "network");
});

await holds("everything else is left to the browser: above all a call to a workstation", async () => {
  const b = browser();
  await b.send("install");
  b.asked.length = 0;
  const left = [
    // the workstation's own server hands the page out at its root, and is called under /v1/ there
    ["http://pages.test/v1/status", {}],
    [SCOPE + "v1/status", {}],
    [SCOPE + "v1/jobs/3/result", {}],
    // a workstation somewhere else, whatever its addresses look like
    ["http://workstation.test:7878/v1/status", {}],
    ["http://workstation.test:7878/helixloom/" + script.slice(2), {}],
    ["http://workstation.test:7878/helixloom/", { mode: "navigate" }],
    // a file of another build, another page of the site, another site
    [SCOPE + "helixloom-0000000000000000.js", {}],
    ["http://pages.test/other/", { mode: "navigate" }],
    ["http://pages.test/" + script.slice(2), {}],
    [SCOPE + "notes/", { mode: "navigate" }],
    ["http://elsewhere.test/helixloom/", { mode: "navigate" }],
    // anything that is not a plain read, even of a kept file
    [urlOf(script), { method: "POST" }],
    [urlOf(script), { method: "HEAD" }],
    [SCOPE, { method: "POST", mode: "navigate" }],
    [urlOf(script) + "?other", {}],
  ];
  for (const [url, how] of left) {
    assert.equal(await b.ask(url, how), "left", url + " " + JSON.stringify(how));
  }
  assert.deepEqual(b.asked, [], "and the worker asked the network for none of them");
});

await holds("when it takes over it removes the caches of this page's older builds, and no other cache", async () => {
  const b = browser();
  const older = "helixloom " + SCOPE + " 0123456789abcdef";
  const others = [
    "another-program",
    // a cache of another HelixLoom page of the same site, in another folder: that page's worker looks after it
    "helixloom http://pages.test/other/ 0123456789abcdef",
    "helixloom-docs",
    "helixloom",
  ];
  [older, ...others].forEach((name) => b.stores.set(name, { held: new Map() }));
  await b.send("install");
  await b.send("activate");
  const names = [...b.stores.keys()];
  assert.ok(!names.includes(older), "the older build's cache is gone");
  for (const other of others) assert.ok(names.includes(other), other + " is left");
  assert.equal(names.length, others.length + 1);
  assert.deepEqual(b.kept().sort(), FILES.map(urlOf).sort(), "and its own is whole");
});

await holds("the manifest is what a browser needs to install the page", async () => {
  const manifest = JSON.parse(readFileSync(manifestPath, "utf8"));
  assert.equal(manifest.name, "HelixLoom");
  assert.equal(manifest.short_name, "HelixLoom");
  assert.equal(manifest.start_url, "./");
  assert.equal(manifest.scope, "./");
  assert.equal(manifest.display, "standalone");
  for (const colour of [manifest.background_color, manifest.theme_color]) assert.match(colour, /^#[0-9a-f]{6}$/i);
  const sizes = manifest.icons.map((icon) => icon.sizes).sort();
  assert.deepEqual(sizes, ["192x192", "512x512"]);
  for (const icon of manifest.icons) {
    assert.equal(icon.type, "image/png");
    assert.ok(!icon.src.startsWith("/") && !icon.src.includes(":"), icon.src + " is beside the manifest");
    assert.ok(FILES.includes("./" + icon.src), icon.src + " is kept by the worker");
  }
});

await holds("the page asks for the worker once it has loaded, and a browser without workers is left in peace", async () => {
  const page = readFileSync(pagePath, "utf8");
  // the page's own script: the one that is not a module
  const scripts = [...page.matchAll(/<script>([\s\S]*?)<\/script>/g)].map((found) => found[1]);
  assert.equal(scripts.length, 1, "one script of the page's own");
  const run = (navigator) => {
    const world = { registered: [], warned: [], onLoad: [] };
    const window = { addEventListener: (type, listener) => type === "load" && world.onLoad.push(listener) };
    const console = { warn: (...said) => world.warned.push(said.join(" ")) };
    vm.runInNewContext(scripts[0], { navigator, window, console });
    return world;
  };
  // a browser with workers: nothing before the page has loaded, the worker beside the page after
  const calls = [];
  const withWorkers = run({ serviceWorker: { register: (address) => (calls.push(address), Promise.resolve()) } });
  assert.deepEqual(calls, []);
  assert.equal(withWorkers.onLoad.length, 1);
  withWorkers.onLoad[0]();
  assert.deepEqual(calls, ["./sw.js"]);
  // a browser that refuses (a private window): said in the console, and nothing breaks
  const refusing = run({ serviceWorker: { register: () => Promise.reject(new Error("refused")) } });
  refusing.onLoad[0]();
  await new Promise((resolve) => setTimeout(resolve, 0));
  assert.equal(refusing.warned.length, 1);
  assert.match(refusing.warned[0], /refused/);
  // a browser without workers, or a page that is not allowed one (plain http to another computer)
  const without = run({});
  assert.equal(without.onLoad.length, 0);
});

console.log(held.map((what) => "held: " + what).join("\n"));
console.log("ALL HELD");
