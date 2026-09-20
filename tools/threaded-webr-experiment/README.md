# Threaded WebR experiment for QCA

This is an isolated feasibility spike. It does not replace DialogQCA's stock
WebR runtime or its published package-library VFS.

The completed benchmark achieved 3.84x speedup with exact output checks.
See [RESULTS.md](RESULTS.md) for measurements and integration limits.

The experiment rebuilds the pinned WebR 0.6.0 runtime with Emscripten pthread
support and a pool of three helper workers. QCA uses four total computation
lanes: the WebR calling worker plus the three helpers. It then builds the current
QCA source as a compatible WebAssembly side module and benchmarks the real
`C_omplexity` entry point with one versus four lanes in Chromium.

The pinned image's optional Cairo libraries are serial WebAssembly archives, so
this spike skips rebuilding that graphics backend while rebuilding the core and
remaining base-package side modules for shared memory. The copied library VFS
can still contain stock serial graphics modules: these are not validated or
suitable for this threaded runtime. Plotting is outside this feasibility boundary.

## Build

Docker must be running and the pinned WebR image must be available.

```sh
bash tools/threaded-webr-experiment/build-runtime.sh
```

All generated runtime and package artifacts remain under the ignored `build/`
directory.

## Run

Point `PLAYWRIGHT_PATH` at an installed Playwright package:

```sh
PLAYWRIGHT_PATH=/absolute/path/to/playwright \
    node tools/threaded-webr-experiment/run.mjs
```

Alternatively, set `PUPPETEER_PATH` to an installed Puppeteer package. The latest
result, including the browser version, is saved to `build/last-result.json`.

The harness serves COOP and COEP headers, requires `crossOriginIsolated`, loads
the QCA side module directly, checks every output against the exact expected
integer result, and reports warmed median timings over six alternating-order
serial/threaded runs. This first benchmark validates the
runtime and package linkage. CCubes and consistent-model DialogQCA workflows
remain the decisive performance benchmark before any application integration.

The experiment lets the WebR worker create native nested browser workers.
WebR 0.6.0 normally replaces `Worker` with an RPC proxy, but that proxy does not
forward Emscripten's pthread module/shared-memory startup message. The native
path keeps that exchange entirely between the WebR worker and its pthread pool.
It also identifies `R.js` as Emscripten's worker entry point, because WebR loads
that script from `webr-worker.js` rather than executing it as the top-level file.
The stock WebR prelude is made safe for that direct worker load by initialising
an empty `Module` object when the surrounding WebR worker has not supplied one.
