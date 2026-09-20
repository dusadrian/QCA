import fs from 'node:fs';

const target = process.argv[2];
if (!target) {
    throw new Error('Usage: node patch-r-runtime.mjs <R.js>');
}

const before = 'var Module = globalThis.Module;';
const after = 'var Module = globalThis.Module || {}; globalThis.Module = Module;';
let source = fs.readFileSync(target, 'utf8');

if (!source.startsWith(before)) {
    throw new Error('Unexpected WebR pre.js header.');
}

source = source.replace(before, after);
fs.writeFileSync(target, source);
