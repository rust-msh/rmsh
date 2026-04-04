// Worker bootstrap: loads the emstudio-worker WASM module and calls the entry point.
import init, { worker_entry } from './emstudio_worker.js';

(async () => {
    await init();
    worker_entry();
})();
