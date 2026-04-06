/**
 * Service Worker for EMStudio Local-First WASM
 * 
 * Features:
 * - Caches all static assets for offline support
 * - Enables Application Shell (AppShell) architecture
 * - Provides version management and cache invalidation
 * - Falls back to network when cache misses
 */

const CACHE_VERSION = 'emstudio-v1';
const CACHE_STATIC = `${CACHE_VERSION}-static`;
const CACHE_DYNAMIC = `${CACHE_VERSION}-dynamic`;

// Assets that must be cached for offline operation
const STATIC_ASSETS = [
  '/',
  '/index.html',
  '/emstudio_main.js',
  '/emstudio_main_bg.wasm',
  '/worker/emstudio_worker.js',
  '/worker/emstudio_worker_bg.wasm',
];

/**
 * Install event: Pre-cache essential assets
 */
self.addEventListener('install', (event) => {
  console.log('[ServiceWorker] Install:', CACHE_STATIC);
  
  event.waitUntil(
    caches
      .open(CACHE_STATIC)
      .then((cache) => {
        console.log('[ServiceWorker] Pre-caching assets...');
        return cache.addAll(STATIC_ASSETS).catch((err) => {
          // Some assets might not exist on first build; that's OK
          console.warn('[ServiceWorker] Pre-cache warning:', err);
        });
      })
      .then(() => {
        // Force new service worker to take control immediately
        return self.skipWaiting();
      })
  );
});

/**
 * Activate event: Clean up old caches
 */
self.addEventListener('activate', (event) => {
  console.log('[ServiceWorker] Activate');
  
  event.waitUntil(
    caches.keys().then((cacheVersions) => {
      return Promise.all(
        cacheVersions.map((version) => {
          if (version !== CACHE_STATIC && version !== CACHE_DYNAMIC) {
            console.log('[ServiceWorker] Deleting old cache:', version);
            return caches.delete(version);
          }
        })
      );
    }).then(() => {
      // Claim all clients immediately
      return self.clients.claim();
    })
  );
});

/**
 * Fetch event: Serve from cache, fall back to network
 * 
 * Strategy: Cache First (for WASM modules) -> Network First (for others)
 */
self.addEventListener('fetch', (event) => {
  const { request } = event;
  const url = new URL(request.url);

  // Skip non-GET requests
  if (request.method !== 'GET') {
    return;
  }

  // Skip non-same-origin requests (CORS)
  if (url.origin !== self.location.origin) {
    return;
  }

  // Determine caching strategy based on request type
  if (request.url.endsWith('.wasm') || request.url.endsWith('.js')) {
    // WASM and JS: Cache FIRST (for fast startup)
    event.respondWith(cacheFirstStrategy(request));
  } else if (request.url.includes('/api/')) {
    // API calls: Network FIRST (for fresh data)
    event.respondWith(networkFirstStrategy(request));
  } else {
    // Everything else: Cache FIRST (static assets)
    event.respondWith(cacheFirstStrategy(request));
  }
});

/**
 * Cache-first strategy: Try cache, fallback to network
 */
async function cacheFirstStrategy(request) {
  try {
    // Try cache first
    const cached = await caches.match(request);
    if (cached) {
      console.log('[ServiceWorker] Cache hit:', request.url);
      return cached;
    }

    // Cache miss, try network
    const response = await fetch(request);

    // Cache successful responses dynamically
    if (response.ok) {
      const cache = await caches.open(CACHE_DYNAMIC);
      cache.put(request, response.clone());
      console.log('[ServiceWorker] Cached from network:', request.url);
    }

    return response;
  } catch (error) {
    console.error('[ServiceWorker] Fetch failed:', error);

    // Return offline page or placeholder
    return new Response('Offline', {
      status: 503,
      statusText: 'Service Unavailable',
    });
  }
}

/**
 * Network-first strategy: Try network, fallback to cache
 */
async function networkFirstStrategy(request) {
  try {
    const response = await fetch(request);
    
    // Cache successful API responses
    if (response.ok) {
      const cache = await caches.open(CACHE_DYNAMIC);
      cache.put(request, response.clone());
    }

    return response;
  } catch (error) {
    console.log('[ServiceWorker] Network failed, trying cache:', request.url);
    
    // Fall back to cache
    const cached = await caches.match(request);
    if (cached) {
      return cached;
    }

    // No cache available
    throw error;
  }
}

/**
 * Handle background sync for offline project saves
 * (Optional: requires Background Sync API support)
 */
self.addEventListener('sync', (event) => {
  if (event.tag === 'sync-projects') {
    console.log('[ServiceWorker] Background sync: projects');
    // Implement sync logic here
  }
});
