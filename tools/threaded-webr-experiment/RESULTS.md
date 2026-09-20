# Browser pthread feasibility result

Verified 2026-09-20 UTC using the completed experimental build.

QCA's real `C_omplexity` C entry point runs in parallel inside the browser.
Median elapsed time fell from **284 ms with one lane to 74 ms with four lanes**:
**3.84x speedup** on this bounded workload.

## Measurement

- Headless Chrome 126.0.6478.182 on macOS; browser reports 12 logical processors.
- WebR 0.6.0 build image, R 4.6.0, Emscripten 5.0.7; rebuilt for shared memory.
- One WebR calling worker plus three native nested pthread workers.
- Local HTTP server with COOP/COEP; `crossOriginIsolated === true`.
- 128 repeated layer-7 complexity calculations for 16 binary conditions.
- Both modes warmed up, then six measurements per mode with alternating order.
- Every returned integer checked against `choose(16, 7) * 2^7 = 1464320`.

| Run | One lane (ms) | Four lanes (ms) |
| --- | ---: | ---: |
| 1 | 292 | 74 |
| 2 | 287 | 74 |
| 3 | 285 | 74 |
| 4 | 281 | 74 |
| 5 | 283 | 73 |
| 6 | 281 | 73 |
| Median | 284 | 74 |

The runner saves raw timing values and browser metadata in the ignored
`build/last-result.json`. These timings exclude startup and module loading.
The one-lane baseline uses the same threaded runtime with QCA restricted to one
lane; it is not a measurement against the deployed stock DialogQCA runtime.

## Checks

The completed Docker build and final browser run both exited successfully.
Native range-dispatch tests passed with pthreads enabled and with the serial
fallback, checking that each item is visited exactly once. The earlier full R
test run in this task passed 2,016 tests with no failures, warnings or skips.
Shell syntax and whitespace checks also passed.

## Boundary and next integration gate

This completes the feasibility spike, not production DialogQCA integration.
The benchmark uses a deliberately parallel-friendly repeated workload. It does
not establish the speedup of truth-table creation, CCubes, minimization, or a
complete analysis, and only one Chromium version has been tested.

Threaded QCA cannot simply be loaded into stock serial WebR. Deployment requires
a compatible runtime and native package/dependency library, including a rebuilt
and tested graphics stack. The experiment skips rebuilding Cairo; copied stock
graphics modules are not suitable evidence of threaded graphics compatibility.
Stock serial Emscripten package builds should explicitly use
`QCA_DISABLE_PTHREAD=1` (or `--disable-pthread`); Emscripten's successful link
probe alone does not establish that the host R runtime supports shared memory.

Before product integration, validate representative DialogQCA workflows and
their full outputs/timings, rebuild dependencies coherently, and test supported
browsers plus selection of the serial runtime where isolation is unavailable.
Only buffer-based native work belongs on helper threads, not concurrent R API
calls. No DialogQCA or shared R-wasm production files were changed or deployed.
