#!/usr/bin/env bash
set -euo pipefail

emcc pocketfft_wasm.cpp \
  -O3 \
  -std=c++17 \
  -s WASM=1 \
  -s MODULARIZE=1 \
  -s EXPORT_ES6=1 \
  -s ENVIRONMENT=web \
  -s ALLOW_MEMORY_GROWTH=1 \
  -s EXPORTED_FUNCTIONS='["_malloc","_free","_dst1_batch_forward","_dst1_batch_inverse"]' \
  -s EXPORTED_RUNTIME_METHODS='["HEAPF64"]' \
  -o pocketfft_wasm.js
