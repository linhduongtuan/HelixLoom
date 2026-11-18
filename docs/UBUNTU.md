# Deploying HelixLoom on Ubuntu

This page takes one Ubuntu computer from nothing to **HelixLoom running as a web app that you open in a browser**, with the
workstation server behind it that does the heavy jobs. The same computer can also run the desktop app (last section). For the
GPU parts (CUDA, PyTorch, folding programs, several cards) this page points to [DEPLOY.md](../DEPLOY.md), which has every option.

**What was tried, and what was not.** HelixLoom has been built, tested and used on a Mac. **Nothing on this page has been run on
Ubuntu**: the commands come from reading the code, the build scripts, the CI file and the libraries' own instructions. They are
what to try first; if a step fails, the message says which, and the notes at the end of each step say what to look at.

```
   your browser (any computer) ──ssh tunnel / VPN / https──▶  Ubuntu computer
                                                              helixloom serve  ──▶  hands out the web page
                                                                                ──▶  computes the jobs it is sent
                                                                                      (alignments, trees, pipelines,
                                                                                       structures, PhyloCNN)
```

You need Ubuntu 22.04 or 24.04 (or later), a user with `sudo`, and about 4 GB of free disk for the build.

## 1. Install what the build needs

```sh
sudo apt update
sudo apt install -y build-essential pkg-config git curl \
    libxcb-render0-dev libxcb-shape0-dev libxcb-xfixes0-dev libxkbcommon-dev libssl-dev
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh      # answer the questions, then open a new shell
rustc --version                                                     # must be a recent stable Rust
```

The server never opens a window, but it is the same program as the desktop app, so the libraries of the window system are
needed to *build* it (the last four are the ones the egui project lists for Linux); a computer with no display is fine to *run*
it on. Use rustup's Rust, not Ubuntu's `rustc` package: the web version needs a Rust target that only rustup installs.

## 2. Get the code and build the program

```sh
git clone https://github.com/linhduongtuan/HelixLoom.git
cd HelixLoom
cargo build --release --no-default-features        # the processor only; the first build takes several minutes
./target/release/helixloom --help 2>&1 | head -5   # it starts
```

`--no-default-features` leaves out the GPU engine (candle), which this computer does not need to serve the web version and to
compute on the processor. **With NVIDIA GPUs** build as in [DEPLOY.md](../DEPLOY.md) step 1.2 (`--features cuda`, or the ordinary
build with PyTorch workers: no CUDA toolkit needed to build).

Optional programs that jobs call when they are installed (a missing one is replaced by the built-in method):

```sh
sudo apt install -y mafft primer3 hmmer
./target/release/helixloom tools            # lists the programs it can call: found or not, version, licence
```

IQ-TREE comes from <https://iqtree.org> or conda ([DEPLOY.md](../DEPLOY.md) step 1.3).

## 3. Build the web version

The web version is a folder of static files (`dist/`) that the server hands out. Either build it here, or build it on your Mac
(where it has been tried) and copy the folder.

**On this Ubuntu computer:**

```sh
rustup target add wasm32-unknown-unknown
grep -A1 'name = "wasm-bindgen"$' Cargo.lock           # the version that Cargo.lock names (0.2.130 at the time of writing)
cargo install --locked wasm-bindgen-cli --version 0.2.130     # use the version printed above; several minutes
scripts/build-web.sh                                    # makes dist/ with a packed copy (.gz) of each file
ls dist
```

The script downloads nothing: when `wasm-bindgen` or the Rust target is missing it says which and stops. **Not tried on Linux**
(it uses `shasum`, or `sha256sum` where only that is there).

**Or on your Mac, then copy:** `scripts/build-web.sh`, then `rsync -a dist/ you@ubuntu-host:HelixLoom/dist/`.

## 4. Try it by hand first

```sh
./target/release/helixloom token > token.txt && chmod 600 token.txt
./target/release/helixloom serve --listen 127.0.0.1:7878 --token-file token.txt --data-dir ./data --web-dir dist
```

It prints where it listens and the name of the computer. From your own computer, open a tunnel and the page:

```sh
ssh -N -L 7878:127.0.0.1:7878 you@ubuntu-host        # leave it running
# then open http://127.0.0.1:7878/ in a browser
```

Go to **Workstation**, paste the token (`cat token.txt` on the Ubuntu computer), press **Connect**. The address is already filled
in. Stop the server with Ctrl+C when it works.

The token is the password: anyone who has it can use the computer. It is never put on a command line (a command line can be read
by other users); it comes from the file.

## 5. Run it as a service that starts at boot

```sh
sudo useradd --system --create-home --home-dir /var/lib/helixloom --shell /usr/sbin/nologin helixloom
sudo mkdir -p /opt/helixloom
sudo cp target/release/helixloom /opt/helixloom/
sudo cp -r dist /opt/helixloom/web                   # only the build goes in this folder: it is handed to anyone who connects
sudo install -d -o helixloom -g helixloom -m 700 /etc/helixloom
sudo -u helixloom sh -c '/opt/helixloom/helixloom token > /etc/helixloom/token' && sudo chmod 600 /etc/helixloom/token
sudo cp deploy/helixloom-server.service /etc/systemd/system/
```

The ready-made unit does not hand out the web page. Add that with a drop-in file, which leaves the unit as it is:

```sh
sudo systemctl edit helixloom-server
```

and in the editor write (the empty `ExecStart=` line clears the unit's own):

```ini
[Service]
ExecStart=
ExecStart=/opt/helixloom/helixloom serve --listen 127.0.0.1:7878 --token-file /etc/helixloom/token --data-dir /var/lib/helixloom/data --gpus auto --web-dir /opt/helixloom/web
```

Then start it:

```sh
sudo systemctl daemon-reload
sudo systemctl enable --now helixloom-server
sudo systemctl status helixloom-server               # active (running)
journalctl -u helixloom-server -f                    # what it says; Ctrl+C leaves
sudo cat /etc/helixloom/token                        # the token to paste into the Workstation page
sudo -u helixloom /opt/helixloom/helixloom status --server 127.0.0.1:7878 --token-file /etc/helixloom/token
```

`--gpus auto` asks `nvidia-smi`; on a computer with no GPU it hands out none and jobs run on the processor. The unit also limits the
server and its jobs to 90 % of the memory, hides `/home`, and runs as its own user; `deploy/helixloom-server.service` explains each line.
The unit was written without an Ubuntu computer to try it on: if systemd complains about a line, `systemctl status` and `journalctl`
say which.

## 6. Reach it from your computers

The server does not encrypt anything, so **never put its port on the open internet**. In the order to prefer:

1. **An SSH tunnel** (step 4): nothing is opened on the Ubuntu computer; the server stays on `127.0.0.1`. The page is at
   `http://127.0.0.1:7878/`, where a browser also lets you **install** it as an app.
2. **A private network** (Tailscale, WireGuard): start the server on that network's address (`--listen 100.64.0.5:7878`, say) and
   open `http://THAT-ADDRESS:7878/`.
3. **Your own trusted network:** `--listen 0.0.0.0:7878` and let only your computers in:
   `sudo ufw allow from 192.168.1.0/24 to any port 7878 proto tcp` (your network's addresses). The page works at
   `http://ubuntu-host:7878/`, but a browser offers to install a page only from `https://` or from the computer itself.
4. **A web address with HTTPS** (a proxy such as Caddy or nginx in front of the server, with a certificate): the only way to
   install the page from another computer, and to use a name such as `https://helix.example.org`. A sketch for Caddy, which gets
   the certificate by itself for a name that points at this computer (**not tried here**):

   ```
   helix.example.org {
       reverse_proxy 127.0.0.1:7878
   }
   ```

   With nginx, `client_max_body_size` must be at least the server's `--max-upload-mb` (512 by default) and `proxy_read_timeout` /
   `proxy_send_timeout` longer than the longest upload or result takes to move ([DEPLOY.md](../DEPLOY.md) section 2.1). Behind a
   public address the token is the only thing between the internet and your computer: prefer 1 or 2 unless you accept that.

The page keeps the address and the token for the **tab** only; *Forget* on the Workstation page removes them
([WEB.md](WEB.md), *The connection, kept for the tab*).

## 7. Update, back up, remove

```sh
cd ~/HelixLoom && git pull
cargo build --release --no-default-features && scripts/build-web.sh
sudo cp target/release/helixloom /opt/helixloom/
sudo rm -rf /opt/helixloom/web && sudo cp -r dist /opt/helixloom/web   # the server reads the folder once, when it starts
sudo systemctl restart helixloom-server
```

Waiting jobs and finished results are kept across a restart; a job that was running at that moment is marked as failed (run it
again). Everything the server keeps is in `/var/lib/helixloom/data` (jobs and results for 7 days, and the answers of pipeline
files in `flow-cache`): copy that folder to back it up. To remove the server:
`sudo systemctl disable --now helixloom-server`, delete `/etc/systemd/system/helixloom-server.service`, `/etc/systemd/system/helixloom-server.service.d`,
`/opt/helixloom`, `/etc/helixloom` and `/var/lib/helixloom`, then `sudo userdel helixloom`.

## 8. The desktop app on an Ubuntu computer with a screen

The packages of step 1 are what the build needs. To *run* it, the computer needs a graphics driver that offers Vulkan (the window
is drawn with wgpu): an NVIDIA computer has it with the driver, otherwise `sudo apt install mesa-vulkan-drivers`. The file
dialogs (File ▸ Open) come from the desktop portal, which an Ubuntu desktop has (`xdg-desktop-portal` with its GNOME or GTK part).

```sh
cargo run --release                          # the app
cargo run --release -- my.gb trace.ab1       # open files
```

Its settings and what it remembers are in `~/.config/helixloom`. **The desktop app has not been tried on Linux**; if the window
does not open, run it from a terminal and read what it says, and send `Help ▸ Copy Diagnostics` with the report if it does open.
How to use it, as an app and in the browser: [USING.md](USING.md).

## 9. If something does not work

| What you see | What to look at |
| --- | --- |
| `cargo build` stops with "could not find system library" | install the `-dev` package it names (`apt search NAME`) |
| `scripts/build-web.sh` says `wasm-bindgen` is missing or the wrong version | install the version `Cargo.lock` names (step 3) |
| the page opens but Workstation says it cannot reach the server | the token is wrong (`sudo cat /etc/helixloom/token`), or the tunnel is not running, or the page was opened from another address than the server's |
| the server refuses to start with `--web-dir` | the folder holds something that is not the build, or the token file or `--data-dir` is inside it (`/opt/helixloom/web` must hold only the build) |
| `systemctl status` says `exited` with a message about the token | the file is missing, not readable by `helixloom`, or shorter than 32 characters: make it with `helixloom token` |
| jobs never start on the GPU | `nvidia-smi` must list the cards; [DEPLOY.md](../DEPLOY.md) steps 1.2 and 1.6 |
