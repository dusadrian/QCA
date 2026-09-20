import { WebR } from './build/runtime/webr.mjs';

const status = document.querySelector('#status');
const report = value => {
    status.textContent = JSON.stringify(value, null, 2);
    globalThis.__qcaThreadedWebRResult = value;
};

try {
    const webR = new WebR({
        baseUrl: new URL('./build/runtime/', window.location.href).href
    });
    status.textContent = 'Initialising WebR';
    await webR.init();
    status.textContent = 'Loading QCA';

    const moduleResponse = await fetch('./build/QCA.so');
    if (!moduleResponse.ok) {
        throw new Error(`QCA.so request failed: ${moduleResponse.status}`);
    }
    await webR.FS.writeFile('/tmp/QCA.so', new Uint8Array(await moduleResponse.arrayBuffer()));

    const result = await webR.evalR(`
        dyn.load("/tmp/QCA.so")
        layers <- rep.int(7L, 128L)
        noflevels <- rep.int(2L, 16L)
        expected <- rep.int(as.integer(choose(16, 7) * 2^7), length(layers))
        run_complexity <- function(threads) {
            Sys.setenv(QCA_NUM_THREADS = as.character(threads))
            timing <- system.time({
                answer <- .Call("C_omplexity", list(16L, layers, noflevels),
                    PACKAGE = "QCA")
            })
            stopifnot(identical(answer, expected))
            unname(timing[["elapsed"]])
        }
        invisible(run_complexity(1L))
        invisible(run_complexity(4L))
        serial <- threaded <- numeric(6L)
        for (i in seq_len(6L)) {
            order <- if (i %% 2L) c(1L, 4L) else c(4L, 1L)
            for (threads in order) {
                elapsed <- run_complexity(threads)
                if (threads == 1L) serial[[i]] <- elapsed else threaded[[i]] <- elapsed
            }
        }
        list(
            serial = serial,
            threaded = threaded,
            identical = TRUE,
            expectedValue = expected[[1L]],
            serialMedian = median(serial),
            threadedMedian = median(threaded),
            speedup = median(serial) / median(threaded)
        )
    `);
    const raw = await result.toJs({ depth: 4 });
    const converted = Object.fromEntries(raw.names.map((name, index) => [
        name, raw.values[index].values.length === 1
            ? raw.values[index].values[0] : raw.values[index].values
    ]));
    report({ ok: converted.identical === true, crossOriginIsolated,
        hardwareConcurrency: navigator.hardwareConcurrency, result: converted });
    webR.close();
} catch (error) {
    report({ ok: false, crossOriginIsolated, error: String(error.stack || error) });
}
