import fs from 'node:fs';

const target = process.argv[2];
if (!target) {
    throw new Error('Usage: node patch-webr-worker.mjs <webr-worker.js>');
}

const before = 'U?.WorkerProxy&&(globalThis.Worker=U.WorkerProxy)';
const after = 'U?.WorkerProxy&&(globalThis.Worker=globalThis.Worker)';
const beforeScript = 'c.locateFile=t=>X.baseUrl+t,';
const afterScript = 'c.mainScriptUrlOrBlob=`${X.baseUrl}R.js`,c.locateFile=t=>X.baseUrl+t,';
let source = fs.readFileSync(target, 'utf8');

for (const marker of [before, beforeScript]) {
    if (source.split(marker).length !== 2) {
        throw new Error(`Expected exactly one WebR worker marker: ${marker}`);
    }
}

source = source.replace(before, after).replace(beforeScript, afterScript);
fs.writeFileSync(target, source);
