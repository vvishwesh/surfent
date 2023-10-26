console.time('Execution Time'); // keeps on tickin'

import { WebR } from 'https://webr.r-wasm.org/latest/webr.mjs'; 

globalThis.webR = new WebR({
    WEBR_URL: "https://webr.r-wasm.org/latest/"
});

console.log("Initializing WebR…")
await globalThis.webR.init();
console.log("WebR Initialized!")

let xfunc = await fetch('scripts/rfuncs.R');
await globalThis.webR.evalR(await xfunc.text());

console.timeEnd('Execution Time');

console.log("WebR Loaded!");


