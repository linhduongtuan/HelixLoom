# The workstation server (`helixloom serve`) with the web version of HelixLoom, for a computer that runs containers.
#
#   docker build -t helixloom-server .
#   helixloom token > token.txt && chmod 600 token.txt           # or any file of 32 characters or more
#   docker run -d --name helixloom -p 127.0.0.1:7878:7878 \
#       -v "$PWD/token.txt:/run/secrets/token:ro" -v helixloom-data:/data helixloom-server
#
# then open http://127.0.0.1:7878/ (from another computer, through an SSH tunnel: DEPLOY.md, section 2). The port is published on
# the loopback address only: the server does not encrypt, so do not publish it to the internet (DEPLOY.md, section 7).
#
# NOT TRIED: this file was written from the build steps of docs/UBUNTU.md and has not been built or run. The server computes on
# the processor; GPUs need the NVIDIA container toolkit and a build made with them (DEPLOY.md, 1.2), which this one is not.

# ---- build the program and the web page
FROM rust:1-bookworm AS build
RUN apt-get update \
    && apt-get install -y --no-install-recommends pkg-config libxcb-render0-dev libxcb-shape0-dev libxcb-xfixes0-dev libxkbcommon-dev libssl-dev \
    && rm -rf /var/lib/apt/lists/*
RUN rustup target add wasm32-unknown-unknown
WORKDIR /src
COPY . .
# the version of wasm-bindgen that Cargo.lock names: the build script of the web page stops on any other
RUN WB="$(grep -A1 'name = "wasm-bindgen"$' Cargo.lock | sed -n 's/^version = "\(.*\)"$/\1/p')" \
    && cargo install --locked wasm-bindgen-cli --version "$WB"
RUN cargo build --release --locked --no-default-features && scripts/build-web.sh

# ---- the server
FROM debian:bookworm-slim
LABEL org.opencontainers.image.source="https://github.com/linhduongtuan/HelixLoom" \
      org.opencontainers.image.description="HelixLoom workstation server and web version" \
      org.opencontainers.image.licenses="MIT"
RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates libxcb1 libxcb-render0 libxcb-shape0 libxcb-xfixes0 libxkbcommon0 libssl3 \
    && rm -rf /var/lib/apt/lists/* \
    && useradd --system --create-home --home-dir /home/helixloom --shell /usr/sbin/nologin helixloom \
    && mkdir -p /data /opt/helixloom \
    && chown helixloom:helixloom /data
COPY --from=build /src/target/release/helixloom /opt/helixloom/helixloom
# only the build of the page goes in this folder: the server hands it to anyone who connects
COPY --from=build /src/dist /opt/helixloom/web
USER helixloom
VOLUME /data
EXPOSE 7878
# it listens on every address of the container; what reaches it is the choice of `docker run -p`
ENTRYPOINT ["/opt/helixloom/helixloom", "serve", "--listen", "0.0.0.0:7878", "--token-file", "/run/secrets/token", "--data-dir", "/data", "--web-dir", "/opt/helixloom/web"]
