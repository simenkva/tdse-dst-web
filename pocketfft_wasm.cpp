#include <emscripten.h>
#include <cmath>

// You provide these headers in ./pocketfft/
#include "pocketfft_hdronly.h"

// pocketfft lives in namespace pocketfft
using pocketfft::shape_t;
using pocketfft::stride_t;

extern "C" {

// Batched in-place DST-I on contiguous rows.
// data: pointer to batch*n doubles laid out row-major
// n: length of each row
// batch: number of rows
EMSCRIPTEN_KEEPALIVE
void dst1_batch_forward(double* data, int n, int batch) {
  // pocketfft expects shape + strides. For 2D array [batch, n] contiguous:
  shape_t shape = { (size_t)batch, (size_t)n };
  stride_t stride = { (ptrdiff_t)(n*sizeof(double)), (ptrdiff_t)sizeof(double) };

  // axes: transform along last axis (axis=1)
  shape_t axes = { 1 };

  // pocketfft dst parameters:
  // dst(in, out, shape, stride_in, stride_out, axes, type, fct, ortho)
  // pocketfft's DST-I ignores ortho; apply orthonormal scaling manually.
  pocketfft::dst(shape, stride, stride, axes, /*type=*/1, data, data, /*fct=*/1.0, /*ortho=*/false);
  const double scale = std::sqrt(2.0 / (double)(n + 1)) * 0.5;
  const size_t total = (size_t)n * (size_t)batch;
  for (size_t i = 0; i < total; i++) data[i] *= scale;
}

EMSCRIPTEN_KEEPALIVE
void dst1_batch_inverse(double* data, int n, int batch) {
  // Orthonormal DST-I is its own inverse.
  dst1_batch_forward(data, n, batch);
}

} // extern "C"
