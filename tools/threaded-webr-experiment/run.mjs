import http from 'node:http';
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

const root = path.dirname(fileURLToPath(import.meta.url));
const timeout = Number(process.env.BROWSER_TIMEOUT || 180000);
const playwrightPath = process.env.PLAYWRIGHT_PATH;
const puppeteerPath = process.env.PUPPETEER_PATH;
let launch;
let newPage;
if (playwrightPath) {
    const { chromium } = await import(pathToFileURL(path.join(playwrightPath, 'index.js')).href);
    launch = () => chromium.launch({ headless: true });
    newPage = browser => browser.newPage();
} else if (puppeteerPath) {
    const { default: puppeteer } = await import(
        pathToFileURL(path.join(puppeteerPath, 'lib/esm/puppeteer/puppeteer.js')).href
    );
    launch = () => puppeteer.launch({ headless: true });
    newPage = browser => browser.newPage();
} else {
    throw new Error('Set PLAYWRIGHT_PATH or PUPPETEER_PATH to an installed browser driver.');
}

const types = new Map([
    ['.html', 'text/html; charset=utf-8'],
    ['.js', 'text/javascript; charset=utf-8'],
    ['.mjs', 'text/javascript; charset=utf-8'],
    ['.wasm', 'application/wasm']
]);
const server = http.createServer((request, response) => {
    const pathname = decodeURIComponent(new URL(request.url, 'http://localhost').pathname);
    const relative = pathname === '/' ? 'benchmark.html' : pathname.slice(1);
    const target = path.resolve(root, relative);
    if (!target.startsWith(`${root}${path.sep}`) || !fs.existsSync(target) || !fs.statSync(target).isFile()) {
        response.writeHead(404).end('Not found');
        return;
    }
    response.writeHead(200, {
        'Content-Type': types.get(path.extname(target)) || 'application/octet-stream',
        'Cross-Origin-Opener-Policy': 'same-origin',
        'Cross-Origin-Embedder-Policy': 'require-corp',
        'Cross-Origin-Resource-Policy': 'same-origin',
        'Cache-Control': 'no-store'
    });
    fs.createReadStream(target).pipe(response);
});

await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
const address = server.address();
const browser = await launch();
try {
    const page = await newPage(browser);
    page.setDefaultTimeout(timeout);
    page.on('console', message => process.stderr.write(`[browser] ${message.text()}\n`));
    page.on('pageerror', error => process.stderr.write(
        `[page error] ${error?.stack || String(error)}\n`
    ));
    page.on('requestfailed', request => process.stderr.write(
        `[request failed] ${request.url()} ${request.failure()?.errorText || ''}\n`
    ));
    page.on('workercreated', worker => process.stderr.write(`[worker] ${worker.url()}\n`));
    await page.goto(`http://127.0.0.1:${address.port}/`, { waitUntil: 'domcontentloaded' });
    if (process.env.WORKER_DEBUG) {
        setTimeout(async() => {
            const workers = typeof page.workers === 'function' ? page.workers() : [];
            const states = await Promise.all(workers.map(async worker => {
                try {
                    return await worker.evaluate(() => ({
                        href: location.href,
                        name: globalThis.name,
                        hasModule: typeof globalThis.Module !== 'undefined',
                        hasPThread: typeof globalThis.PThread !== 'undefined',
                        hasOnMessage: typeof globalThis.onmessage,
                        hasWasmMemory: typeof globalThis.wasmMemory !== 'undefined',
                        pthreadUnused: globalThis.PThread?.unusedWorkers?.length,
                        pthreadLoaded: globalThis.PThread?.unusedWorkers?.map(worker => ({
                            loaded: worker.loaded,
                            uuid: worker.uuid
                        }))
                    }));
                } catch (error) {
                    return { href: worker.url(), error: String(error) };
                }
            }));
            process.stderr.write(`[worker states] ${JSON.stringify(states)}\n`);
        }, 10000);
    }
    try {
        await page.waitForFunction(() => globalThis.__qcaThreadedWebRResult);
    } catch (error) {
        const state = await page.evaluate(() => ({
            crossOriginIsolated,
            status: document.querySelector('#status')?.textContent,
            worker: typeof Worker,
            hardwareConcurrency: navigator.hardwareConcurrency
        }));
        process.stderr.write(`[timeout state] ${JSON.stringify(state)}\n`);
        throw error;
    }
    const result = await page.evaluate(() => globalThis.__qcaThreadedWebRResult);
    result.browserVersion = await browser.version();
    result.measuredAt = new Date().toISOString();
    fs.writeFileSync(path.join(root, 'build', 'last-result.json'), `${JSON.stringify(result, null, 2)}\n`);
    process.stdout.write(`${JSON.stringify(result, null, 2)}\n`);
    if (!result.ok || !result.crossOriginIsolated) process.exitCode = 1;
} finally {
    await browser.close();
    server.close();
}
