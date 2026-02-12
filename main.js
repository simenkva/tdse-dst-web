// 2D TDSE split-operator using DST-I (Dirichlet) via pocketfft WASM.
// Units: ħ = 1, mass = 1.
// H = -1/2 ∇^2 + V(x,y)
// Preset potentials loaded from grayscale PNGs and edited with a brush.
//
// Grid: N interior points per axis, boundaries at x=±L, y=±L with ψ=0 there.
// dx = 2L/(N+1); x_i = -L + (i+1)dx, i=0..N-1.
//
// Kinetic eigenvalues for DST-I modes m=1..N: k_m = mπ/(2L),
// T_mn = 0.5*(k_m^2 + k_n^2).

let wasm = null;

async function loadWasm() {
  // Emscripten MODULARIZE=1 build exports a factory function.
  const modFactory = (await import("./pocketfft_wasm.js")).default;
  wasm = await modFactory();
  return wasm;
}

// ---------- UI visibility config ----------
const UI_CONFIG = {
  run: true,
  grid: true,
  domain: true,
  timestep: true,
  steps: false,
  preset: true,
  symbolicPotential: true,
  integrator: false,
  visualization: true,
  brushSize: true,
  brushHardness: false,
  brushValue: true,
  x0: true,
  y0: true,
  vMag: true,
  vAngle: true,
  sigma: true,
  mask: false,
  renorm: false,
};

// ---------- Default values ----------
const DEFAULTS = {
  k: 9,
  L: 10,
  dt: 0.0025,
  steps: 1,
  presetIndex: 0,
  integrator: "lie", // "strang" or "lie"
  visMode: "phase",
  brushSize: 1.0,
  brushHardness: 0.7,
  brushValue: 100,
  x0: -5.0,
  y0: 0,
  vMag: 5,
  vAngle: 0,
  sigma: 1,
  mask: 0.6,
  renorm: false,
  symbolicExpr: "",
  symbolicExampleIndex: 0,
};

// Slider ranges (min/max/step) for inputs that should be configurable.
const DEFAULT_RANGES = {
  k: { min: 6, max: 9, step: 1 },
  L: { min: 5, max: 25, step: 0.1 },
  dt: { min: 0.0001, max: 0.01, step: 0.0001 },
  steps: { min: 1, max: 10, step: 1 },
  brushSize: { min: 0.1, max: 6, step: 0.05 },
  brushHardness: { min: 0, max: 1, step: 0.01 },
  brushValue: { min: -5, max: 5, step: 0.05 },
  x0: { min: -6, max: 6, step: 0.05 },
  y0: { min: -6, max: 6, step: 0.05 },
  vMag: { min: 0, max: 20, step: 0.1 },
  vAngle: { min: 0, max: 360, step: 1 },
  sigma: { min: 0.2, max: 4, step: 0.05 },
  mask: { min: 0, max: 3, step: 0.05 },
};

// ---------- UI ----------
const el = {
  canvas: document.getElementById("canvas"),
  btnToggle: document.getElementById("btnToggle"),
  btnResetPsi: document.getElementById("btnResetPsi"),
  btnResetPotential: document.getElementById("btnResetPotential"),
  btnResetAll: document.getElementById("btnResetAll"),
  btnIntegrator: document.getElementById("btnIntegrator"),
  k: document.getElementById("k"),
  L: document.getElementById("L"),
  dt: document.getElementById("dt"),
  steps: document.getElementById("steps"),
  preset: document.getElementById("preset"),
  symbolicExprExample: document.getElementById("symbolicExprExample"),
  symbolicExpr: document.getElementById("symbolicExpr"),
  btnEvalSymbolic: document.getElementById("btnEvalSymbolic"),
  btnSaveSymbolic: document.getElementById("btnSaveSymbolic"),
  btnDeleteSymbolic: document.getElementById("btnDeleteSymbolic"),
  visMode: document.getElementById("visMode"),
  mask: document.getElementById("mask"),
  brushSize: document.getElementById("brushSize"),
  brushHardness: document.getElementById("brushHardness"),
  brushValue: document.getElementById("brushValue"),
  x0: document.getElementById("x0"),
  y0: document.getElementById("y0"),
  vMag: document.getElementById("vMag"),
  vAngle: document.getElementById("vAngle"),
  sigma: document.getElementById("sigma"),
  kVal: document.getElementById("kVal"),
  LVal: document.getElementById("LVal"),
  dtVal: document.getElementById("dtVal"),
  stepsVal: document.getElementById("stepsVal"),
  presetVal: document.getElementById("presetVal"),
  integratorVal: document.getElementById("integratorVal"),
  presetInfo: document.getElementById("presetInfo"),
  maskVal: document.getElementById("maskVal"),
  brushSizeVal: document.getElementById("brushSizeVal"),
  brushHardnessVal: document.getElementById("brushHardnessVal"),
  brushValueVal: document.getElementById("brushValueVal"),
  x0Val: document.getElementById("x0Val"),
  y0Val: document.getElementById("y0Val"),
  vMagVal: document.getElementById("vMagVal"),
  vAngleVal: document.getElementById("vAngleVal"),
  sigmaVal: document.getElementById("sigmaVal"),
  renorm: document.getElementById("renorm"),
  stats: document.getElementById("stats"),
};

function applyDefaultsToElements() {
  if (el.k) { el.k.value = DEFAULTS.k; el.k.defaultValue = String(DEFAULTS.k); }
  if (el.L) { el.L.value = DEFAULTS.L; el.L.defaultValue = String(DEFAULTS.L); }
  if (el.dt) { el.dt.value = DEFAULTS.dt; el.dt.defaultValue = String(DEFAULTS.dt); }
  if (el.steps) { el.steps.value = DEFAULTS.steps; el.steps.defaultValue = String(DEFAULTS.steps); }
  if (el.visMode && DEFAULTS.visMode) {
    for (const opt of el.visMode.options) opt.defaultSelected = false;
    const idx = Array.from(el.visMode.options).findIndex((o) => o.value === DEFAULTS.visMode);
    if (idx >= 0) el.visMode.options[idx].defaultSelected = true;
    el.visMode.value = DEFAULTS.visMode;
  }
  if (el.brushSize) { el.brushSize.value = DEFAULTS.brushSize; el.brushSize.defaultValue = String(DEFAULTS.brushSize); }
  if (el.brushHardness) { el.brushHardness.value = DEFAULTS.brushHardness; el.brushHardness.defaultValue = String(DEFAULTS.brushHardness); }
  if (el.brushValue) { el.brushValue.value = DEFAULTS.brushValue; el.brushValue.defaultValue = String(DEFAULTS.brushValue); }
  if (el.x0) { el.x0.value = DEFAULTS.x0; el.x0.defaultValue = String(DEFAULTS.x0); }
  if (el.y0) { el.y0.value = DEFAULTS.y0; el.y0.defaultValue = String(DEFAULTS.y0); }
  if (el.vMag) { el.vMag.value = DEFAULTS.vMag; el.vMag.defaultValue = String(DEFAULTS.vMag); }
  if (el.vAngle) { el.vAngle.value = DEFAULTS.vAngle; el.vAngle.defaultValue = String(DEFAULTS.vAngle); }
  if (el.sigma) { el.sigma.value = DEFAULTS.sigma; el.sigma.defaultValue = String(DEFAULTS.sigma); }
  if (el.mask) { el.mask.value = DEFAULTS.mask; el.mask.defaultValue = String(DEFAULTS.mask); }
  if (el.symbolicExpr) {
    el.symbolicExpr.value = DEFAULTS.symbolicExpr || "";
    el.symbolicExpr.defaultValue = DEFAULTS.symbolicExpr || "";
  }
  if (el.renorm) {
    el.renorm.checked = !!DEFAULTS.renorm;
    el.renorm.defaultChecked = !!DEFAULTS.renorm;
  }
}

function applyRangesToElements() {
  for (const [key, cfg] of Object.entries(DEFAULT_RANGES)) {
    const ctrl = el[key];
    if (!ctrl || ctrl.tagName !== "INPUT" || ctrl.type !== "range") continue;
    if (cfg.min !== undefined) ctrl.min = String(cfg.min);
    if (cfg.max !== undefined) ctrl.max = String(cfg.max);
    if (cfg.step !== undefined) ctrl.step = String(cfg.step);
  }
}

applyRangesToElements();
applyDefaultsToElements();

function updateLabels(state) {
  el.kVal.textContent = `${state.k} (N=${state.N})`;
  el.LVal.textContent = `${state.L.toFixed(2)}`;
  el.dtVal.textContent = `${state.dt.toFixed(4)}`;
  el.stepsVal.textContent = `${state.stepsPerFrame}`;
  const potentialSource = state.symbolicPotentialActive
    ? "Symbolic expression"
    : (state.preset ? state.preset.name : "");
  el.presetVal.textContent = potentialSource;
  if (el.integratorVal) el.integratorVal.textContent = getIntegratorLabel(state.integrator);
  el.maskVal.textContent = `${state.maskExtent.toFixed(2)}`;
  el.brushSizeVal.textContent = `${state.brushSize.toFixed(2)}`;
  el.brushHardnessVal.textContent = `${state.brushHardness.toFixed(2)}`;
  el.brushValueVal.textContent = `${state.brushValue.toFixed(2)}`;
  el.x0Val.textContent = `${state.x0.toFixed(2)}`;
  el.y0Val.textContent = `${state.y0.toFixed(2)}`;
  el.vMagVal.textContent = `${state.vMag.toFixed(2)}`;
  el.vAngleVal.textContent = `${state.vAngleDeg.toFixed(0)}`;
  el.sigmaVal.textContent = `${state.sigma.toFixed(2)}`;
}

function clamp(x, a, b) { return Math.max(a, Math.min(b, x)); }

function getIntegratorLabel(mode) {
  if (mode === "lie") return "Lie (T - V)";
  return "Strang (V/2 - T - V/2)";
}

const POTENTIAL_PRESETS = [
  // Update these entries to match PNGs in ./potentials (grayscale: 0 => Vmin, 255 => Vmax).
  { name: "Flat", file: null, vMin: 0, vMax: 100 },
  { name: "Double slit experiment", file: "potentials/double_slit.png", vMin: 0, vMax: 100 },
  { name: "Cylinder scattering", file: "potentials/cylinder.png", vMin: 0, vMax: 100 },
  { name: "Stadium", file: "potentials/stadium.png", vMin: 0, vMax: 100 },
  { name: "Crystal scattering 1", file: "potentials/crystal.png", vMin: 0, vMax: 100 },
  { name: "Crystal scattering 2", file: "potentials/crystal2.png", vMin: 0, vMax: 100 },
  { name: "Crystal scattering 3", file: "potentials/crystal3.png", vMin: 0, vMax: 100 },
  { name: "Crystal scattering 4", file: "potentials/crystal4.png", vMin: 0, vMax: 100 },
  { name: "Waveguide", file: "potentials/guides.png", vMin: 0, vMax: 100 },
  { name: "Mario", file: "potentials/mario.png", vMin: 0, vMax: 50 },
  { name: "Gravity", file: "potentials/gravity.png", vMin: 0, vMax: 500 },
  { name: "Harmonic oscillator", file: "potentials/harmonic_osc.png", vMin: 0, vMax: 600 },
];

// Configurable examples for symbolic potential expressions.
const SYMBOLIC_POTENTIAL_EXAMPLES = [
  { name: "Harmonic bowl", expression: "0.5 * (x*x + y*y)" },
  { name: "Ring barrier", expression: "120 * exp(-((r - 4.2)*(r - 4.2)) / 0.35)" },
  { name: "Double slit wall", expression: "abs(x) < 0.22 && (abs(y) > 1.1 || abs(y) < 0.35) ? 160 : 0" },
  { name: "parabola", expression: "100 * (x - 9 > -.5 * y**2)" },
  { name: "Tilted plane", expression: "35 + 7*x - 4*y" },
  { name: "Periodic lattice", expression: "30 * (sin(2.6*x)*sin(2.6*x) + sin(2.6*y)*sin(2.6*y))" },
];
const SYMBOLIC_SAVED_STORAGE_KEY = "tdse2d.symbolicExamples.v1";

function loadSavedSymbolicExamples() {
  try {
    const raw = localStorage.getItem(SYMBOLIC_SAVED_STORAGE_KEY);
    if (!raw) return [];
    const parsed = JSON.parse(raw);
    if (!Array.isArray(parsed)) return [];
    const normalized = [];
    for (const item of parsed) {
      if (!item || typeof item !== "object") continue;
      const name = typeof item.name === "string" ? item.name.trim() : "";
      const expression = typeof item.expression === "string" ? item.expression.trim() : "";
      if (!name || !expression) continue;
      normalized.push({ name, expression });
    }
    return normalized;
  } catch {
    return [];
  }
}

function saveSymbolicExamplesToStorage(examples) {
  try {
    localStorage.setItem(SYMBOLIC_SAVED_STORAGE_KEY, JSON.stringify(examples));
    return true;
  } catch {
    return false;
  }
}

async function loadPresetImage(preset) {
  if (!preset.file) return null;
  if (preset.image) return preset.image;
  const img = new Image();
  img.src = preset.file;
  await new Promise((resolve, reject) => {
    img.onload = resolve;
    img.onerror = reject;
  });
  preset.image = img;
  return img;
}

function buildSymbolicPotentialEvaluator(expression) {
  const raw = String(expression ?? "").trim();
  if (!raw) throw new Error("Expression is empty.");

  // Treat '^' as exponent to support common symbolic syntax.
  const jsExpr = raw.replace(/\^/g, "**");

  const fn = new Function(
    "x", "y", "r", "theta", "L", "N", "dx", "i", "j", "pi", "e", "clamp",
    `"use strict";
const { abs, acos, acosh, asin, asinh, atan, atan2, atanh, cbrt, ceil, cos, cosh, exp, expm1, floor, hypot, log, log10, log1p, log2, max, min, pow, round, sign, sin, sinh, sqrt, tan, tanh, trunc } = Math;
return (${jsExpr});`
  );
  return { raw, fn };
}

// ---------- Color mapping (HSV: hue=phase, value=|psi|) ----------
function hsvToRgb(h, s, v) {
  // h in [0,1)
  const i = Math.floor(h * 6);
  const f = h * 6 - i;
  const p = v * (1 - s);
  const q = v * (1 - f * s);
  const t = v * (1 - (1 - f) * s);
  let r, g, b;
  switch (i % 6) {
    case 0: r=v; g=t; b=p; break;
    case 1: r=q; g=v; b=p; break;
    case 2: r=p; g=v; b=t; break;
    case 3: r=p; g=q; b=v; break;
    case 4: r=t; g=p; b=v; break;
    case 5: r=v; g=p; b=q; break;
  }
  return [r, g, b];
}

function lerp(a, b, t) { return a + (b - a) * t; }

// ---------- Simulation state ----------
class TDSE2D {
  constructor() {
    this.running = false;

    // defaults
    this.k = parseInt(el.k.value, 10);
    this.N = (1 << this.k) - 1;
    this.L = parseFloat(el.L.value);
    this.dt = parseFloat(el.dt.value);
    this.stepsPerFrame = parseInt(el.steps.value, 10);
    this.renormEnabled = !!el.renorm.checked;
    this.maskExtent = parseFloat(el.mask.value);
    this.integrator = DEFAULTS.integrator || "strang";

    // preset + brush
    this.presetIndex = 0;
    this.preset = POTENTIAL_PRESETS[this.presetIndex];
    this.visMode = el.visMode ? el.visMode.value : "phase";
    this.brushSize = parseFloat(el.brushSize.value);
    this.brushHardness = parseFloat(el.brushHardness.value);
    this.brushValue = parseFloat(el.brushValue.value);
    this.symbolicExpression = el.symbolicExpr ? el.symbolicExpr.value.trim() : "";
    this.symbolicPotentialActive = false;
    this.symbolicRangeMin = 0;
    this.symbolicRangeMax = 0;

    // gaussian initial condition
    this.x0 = parseFloat(el.x0.value);
    this.y0 = parseFloat(el.y0.value);
    this.vMag = parseFloat(el.vMag.value);
    this.vAngleDeg = parseFloat(el.vAngle.value);
    this.sigma = parseFloat(el.sigma.value);

    // drawing state
    this.painting = false;

    // arrays
    this.dx = 0;
    this.x = null;
    this.y = null;

    this.psiRe = null;
    this.psiIm = null;

    this.VBase = null; // base potential from preset
    this.V = null;     // editable potential V(x,y)

    this.kx2 = null;  // kx^2 for modes m=1..N
    this.ky2 = null;
    this.maskX = null;
    this.maskY = null;

    this.KRe = null;  // kinetic phase exp(-i dt * 0.5*(kx^2+ky^2))
    this.KIm = null;  // stored per (m,n), size N*N

    // transpose buffers
    this.tmpRe = null;
    this.tmpIm = null;

    // wasm buffers
    this.ptrRe = 0;
    this.ptrIm = 0;
    this.nBytes = 0;

    this.t = 0;
    this.frameCount = 0;
    this.lastFpsTime = performance.now();
    this.fps = 0;
    this.loadingPreset = false;
    this.presetToken = 0;
    this.overlayUntil = 0;
    this.hoverActive = false;
    this.hoverX = 0;
    this.hoverY = 0;
  }

  resizeCanvas() {
    const dpr = window.devicePixelRatio || 1;
    const rect = el.canvas.getBoundingClientRect();
    el.canvas.width = Math.floor(rect.width * dpr);
    el.canvas.height = Math.floor(rect.height * dpr);
  }

  async allocate() {
    this.N = (1 << this.k) - 1;

    const N = this.N;
    const NN = N * N;

    this.dx = (2 * this.L) / (N + 1);
    this.x = new Float64Array(N);
    this.y = new Float64Array(N);
    for (let i = 0; i < N; i++) {
      this.x[i] = -this.L + (i + 1) * this.dx;
      this.y[i] = -this.L + (i + 1) * this.dx;
    }

    this.psiRe = new Float64Array(NN);
    this.psiIm = new Float64Array(NN);

    this.VBase = new Float64Array(NN);
    this.V = new Float64Array(NN);

    this.tmpRe = new Float64Array(NN);
    this.tmpIm = new Float64Array(NN);

    // Precompute k^2 for DST-I modes m=1..N: k_m = m*pi/(2L)
    const factor = Math.PI / (2 * this.L);
    this.kx2 = new Float64Array(N);
    this.ky2 = new Float64Array(N);
    for (let m = 1; m <= N; m++) {
      const km = m * factor;
      this.kx2[m - 1] = km * km;
      this.ky2[m - 1] = km * km;
    }

    // Allocate kinetic phase arrays
    this.KRe = new Float64Array(NN);
    this.KIm = new Float64Array(NN);

    // Allocate wasm memory for re/im arrays (double)
    this.nBytes = NN * 8;
    if (this.ptrRe) wasm._free(this.ptrRe);
    if (this.ptrIm) wasm._free(this.ptrIm);
    this.ptrRe = wasm._malloc(this.nBytes);
    this.ptrIm = wasm._malloc(this.nBytes);

    // Fill V and initial psi
    await this.applyPresetToGrid();
    if (this.symbolicPotentialActive && this.symbolicExpression) {
      const applied = this.applySymbolicExpression(this.symbolicExpression, { updateInput: false });
      if (!applied.ok) {
        this.symbolicPotentialActive = false;
      }
    }
    this.initPsiGaussian();
    this.buildKineticPhase(); // depends on dt
    this.buildMask();
    this.dstRoundTripGain = measureDstRoundTripGain(this.N);
  }

  async applyPresetToGrid(token = null) {
    const preset = this.preset;
    const N = this.N;
    const NN = N * N;
    const vMin = preset ? preset.vMin : 0;
    const vMax = preset ? preset.vMax : 0;

    this.loadingPreset = true;
    el.presetInfo.textContent = preset
      ? `Vmin=${vMin.toFixed(2)}, Vmax=${vMax.toFixed(2)}`
      : "";

    let img = null;
    try {
      img = preset ? await loadPresetImage(preset) : null;
    } catch {
      img = null;
    }

    if (token !== null && token !== this.presetToken) {
      this.loadingPreset = false;
      return;
    }

    if (!img) {
      if (preset && preset.file) {
        el.presetInfo.textContent =
          `Missing ${preset.file}. Using flat V=${vMin.toFixed(2)}.`;
      }
      this.VBase.fill(vMin);
    } else {
      const off = document.createElement("canvas");
      off.width = N;
      off.height = N;
      const octx = off.getContext("2d", { willReadFrequently: true });
      octx.drawImage(img, 0, 0, N, N);
      const data = octx.getImageData(0, 0, N, N).data;
      const scale = (vMax - vMin) / 255;
      for (let j = 0; j < N; j++) {
        for (let i = 0; i < N; i++) {
          const imgIdx = i + N * j;
          const dstIdx = i + N * (N - 1 - j);
          const gray = data[4 * imgIdx];
          this.VBase[dstIdx] = vMin + gray * scale;
        }
      }
    }

    this.V.set(this.VBase);
    this.loadingPreset = false;
  }

  initPsiGaussian() {
    // psi = C * exp(-((x-x0)^2+(y-y0)^2)/(2*sigma^2)) * exp(i (v·(x-x0,y-y0)))
    if (!this.psiRe || !this.psiIm || !this.x || !this.y) return;
    const N = this.N;
    const theta = (this.vAngleDeg * Math.PI) / 180;
    const vx = this.vMag * Math.cos(theta);
    const vy = this.vMag * Math.sin(theta);
    const sig2 = this.sigma * this.sigma;
    const inv2sig2 = sig2 > 0 ? 1 / (2 * sig2) : 1;
    for (let j = 0; j < N; j++) {
      const y = this.y[j];
      const dy = y - this.y0;
      for (let i = 0; i < N; i++) {
        const x = this.x[i];
        const dx = x - this.x0;
        const idx = i + N * j;
        const amp = Math.exp(-(dx * dx + dy * dy) * inv2sig2);
        const phase = vx * dx + vy * dy;
        this.psiRe[idx] = amp * Math.cos(phase);
        this.psiIm[idx] = amp * Math.sin(phase);
      }
    }
    this.t = 0;
  }

  buildKineticPhase() {
    // K(m,n) = exp(-i dt * 0.5*(kx^2 + ky^2))
    const N = this.N;
    for (let n = 0; n < N; n++) {
      const ky2 = this.ky2[n];
      for (let m = 0; m < N; m++) {
        const kx2 = this.kx2[m];
        const E = 0.5 * (kx2 + ky2);
        const theta = -this.dt * E; // ħ=1
        const idx = m + N * n;
        this.KRe[idx] = Math.cos(theta);
        this.KIm[idx] = Math.sin(theta);
      }
    }
  }

  async updateParamsFromUI() {
    const newK = parseInt(el.k.value, 10);
    const newN = (1 << newK) - 1;
    const newL = parseFloat(el.L.value);
    const newDt = parseFloat(el.dt.value);
    const newSteps = parseInt(el.steps.value, 10);
    const newRenorm = !!el.renorm.checked;

    const needRealloc = (newK !== this.k) || (newL !== this.L);

    this.k = newK;
    this.N = newN;
    this.L = newL;
    this.dt = newDt;
    this.stepsPerFrame = Number.isFinite(newSteps) ? newSteps : 1;
    const newMask = parseFloat(el.mask.value);
    this.renormEnabled = newRenorm;
    this.maskExtent = newMask;

    this.syncRangesForL();
    updateLabels(this);

    if (needRealloc) {
      await this.allocate();
    } else {
      // dt changed => kinetic phase must be rebuilt
      this.buildKineticPhase();
      this.buildMask();
    }
  }

  syncRangesForL() {
    const N = (1 << this.k) - 1;
    const L = this.L;
    const dx = (2 * L) / (N + 1);

    const minBrush = Math.max(dx, 0.02);
    const maxBrush = 2 * L;
    el.brushSize.min = minBrush.toFixed(3);
    el.brushSize.max = maxBrush.toFixed(2);
    this.brushSize = clamp(parseFloat(el.brushSize.value), minBrush, maxBrush);
    el.brushSize.value = this.brushSize;

    el.x0.min = (-L).toFixed(2);
    el.x0.max = (L).toFixed(2);
    el.y0.min = (-L).toFixed(2);
    el.y0.max = (L).toFixed(2);
    this.x0 = clamp(parseFloat(el.x0.value), -L, L);
    this.y0 = clamp(parseFloat(el.y0.value), -L, L);
    el.x0.value = this.x0;
    el.y0.value = this.y0;

    const minSigma = 0.2;
    const maxSigma = Math.max(minSigma, L);
    el.sigma.min = minSigma.toFixed(2);
    el.sigma.max = maxSigma.toFixed(2);
    this.sigma = clamp(parseFloat(el.sigma.value), minSigma, maxSigma);
    el.sigma.value = this.sigma;

    const maxMask = Math.max(0, L);
    el.mask.min = "0";
    el.mask.max = maxMask.toFixed(2);
    this.maskExtent = clamp(parseFloat(el.mask.value), 0, maxMask);
    el.mask.value = this.maskExtent;
  }

  getCurrentPotentialRange() {
    if (this.symbolicPotentialActive && Number.isFinite(this.symbolicRangeMin) && Number.isFinite(this.symbolicRangeMax)) {
      return { min: this.symbolicRangeMin, max: this.symbolicRangeMax };
    }

    const presetMin = this.preset ? this.preset.vMin : 0;
    const presetMax = this.preset ? this.preset.vMax : 0;
    if (Number.isFinite(presetMin) && Number.isFinite(presetMax) && presetMax > presetMin) {
      return { min: presetMin, max: presetMax };
    }

    const src = this.VBase && this.VBase.length ? this.VBase : this.V;
    if (!src || !src.length) return { min: 0, max: 0 };

    let min = Infinity;
    let max = -Infinity;
    for (let idx = 0; idx < src.length; idx++) {
      const v = src[idx];
      if (!Number.isFinite(v)) continue;
      if (v < min) min = v;
      if (v > max) max = v;
    }
    if (!Number.isFinite(min) || !Number.isFinite(max)) return { min: 0, max: 0 };
    return { min, max };
  }

  syncBrushValueRange() {
    const { min: rangeMin, max: rangeMax } = this.getCurrentPotentialRange();
    let min = rangeMin;
    let max = rangeMax;
    if (Math.abs(rangeMax - rangeMin) < 1e-6) {
      min = rangeMin - 5;
      max = rangeMax + 5;
    }
    el.brushValue.min = min.toFixed(2);
    el.brushValue.max = max.toFixed(2);
    this.brushValue = clamp(parseFloat(el.brushValue.value), min, max);
    el.brushValue.value = this.brushValue;
  }

  applySymbolicExpression(expression, options = {}) {
    const { updateInput = true } = options;
    if (!this.VBase || !this.V || !this.x || !this.y) {
      el.presetInfo.textContent = "Grid is not initialized yet.";
      return { ok: false, error: "Grid is not initialized yet." };
    }

    let compiled = null;
    try {
      compiled = buildSymbolicPotentialEvaluator(expression);
    } catch (err) {
      const message = `Expression parse error: ${err.message}`;
      el.presetInfo.textContent = message;
      return { ok: false, error: message };
    }

    const N = this.N;
    const next = new Float64Array(N * N);
    let min = Infinity;
    let max = -Infinity;
    try {
      for (let j = 0; j < N; j++) {
        const y = this.y[j];
        for (let i = 0; i < N; i++) {
          const x = this.x[i];
          const r = Math.hypot(x, y);
          const theta = Math.atan2(y, x);
          const idx = i + N * j;
          const value = Number(compiled.fn(x, y, r, theta, this.L, N, this.dx, i, j, Math.PI, Math.E, clamp));
          if (!Number.isFinite(value)) {
            throw new Error(`Non-finite value at i=${i}, j=${j}`);
          }
          next[idx] = value;
          if (value < min) min = value;
          if (value > max) max = value;
        }
      }
    } catch (err) {
      const message = `Expression evaluation error: ${err.message}`;
      el.presetInfo.textContent = message;
      return { ok: false, error: message };
    }

    this.VBase.set(next);
    this.V.set(next);
    this.symbolicExpression = compiled.raw;
    this.symbolicPotentialActive = true;
    this.symbolicRangeMin = min;
    this.symbolicRangeMax = max;
    if (updateInput && el.symbolicExpr) {
      el.symbolicExpr.value = this.symbolicExpression;
    }
    this.syncBrushValueRange();
    updateLabels(this);
    el.presetInfo.textContent = `Symbolic potential active (Vmin=${min.toFixed(2)}, Vmax=${max.toFixed(2)}).`;
    return { ok: true };
  }

  updateBrushFromUI() {
    this.brushSize = parseFloat(el.brushSize.value);
    this.brushHardness = parseFloat(el.brushHardness.value);
    this.brushValue = parseFloat(el.brushValue.value);
    updateLabels(this);
  }

  updateGaussianFromUI() {
    this.x0 = parseFloat(el.x0.value);
    this.y0 = parseFloat(el.y0.value);
    this.vMag = parseFloat(el.vMag.value);
    this.vAngleDeg = parseFloat(el.vAngle.value);
    this.sigma = parseFloat(el.sigma.value);
    updateLabels(this);
    this.initPsiGaussian();
    this.overlayUntil = performance.now() + 700;
  }

  buildMask() {
    const N = this.N;
    const L = this.L;
    const extent = clamp(this.maskExtent, 0, L);
    this.maskX = new Float64Array(N);
    this.maskY = new Float64Array(N);
    if (extent <= 0) {
      this.maskX.fill(1);
      this.maskY.fill(1);
      return;
    }
    const invExtent = 1 / extent;
    for (let i = 0; i < N; i++) {
      const d = L - Math.abs(this.x[i]);
      const t = clamp(d * invExtent, 0, 1);
      const w = Math.sin(0.5 * Math.PI * t);
      this.maskX[i] = w;
      this.maskY[i] = w;
    }
  }

  async setPreset(index) {
    const next = clamp(index, 0, POTENTIAL_PRESETS.length - 1);
    this.presetIndex = next;
    this.preset = POTENTIAL_PRESETS[next];
    this.symbolicPotentialActive = false;
    const token = ++this.presetToken;
    el.preset.value = String(next);
    this.syncBrushValueRange();
    updateLabels(this);
    if (this.VBase && this.V) {
      await this.applyPresetToGrid(token);
    }
  }

  resetPotential() {
    if (this.V && this.VBase) this.V.set(this.VBase);
  }

  paintAt(x0, y0) {
    if (!this.V || !this.x || !this.y) return;
    const N = this.N;
    const R = this.brushSize * 0.5;
    if (R <= 0) return;
    const dx = this.dx;
    const xMin = x0 - R;
    const xMax = x0 + R;
    const yMin = y0 - R;
    const yMax = y0 + R;
    const iMin = clamp(Math.floor((xMin + this.L) / dx) - 1, 0, N - 1);
    const iMax = clamp(Math.ceil((xMax + this.L) / dx) - 1, 0, N - 1);
    const jMin = clamp(Math.floor((yMin + this.L) / dx) - 1, 0, N - 1);
    const jMax = clamp(Math.ceil((yMax + this.L) / dx) - 1, 0, N - 1);

    const inner = R * clamp(this.brushHardness, 0, 1);
    const denom = Math.max(R - inner, 1e-6);
    for (let j = jMin; j <= jMax; j++) {
      const dy = this.y[j] - y0;
      for (let i = iMin; i <= iMax; i++) {
        const dxw = this.x[i] - x0;
        const r = Math.hypot(dxw, dy);
        if (r > R) continue;
        let w = 1;
        if (r > inner) {
          const t = (r - inner) / denom;
          w = 1 - t * t * (3 - 2 * t);
        }
        const idx = i + N * j;
        const v = this.V[idx];
        this.V[idx] = v + (this.brushValue - v) * w;
      }
    }
  }

  setIntegrator(mode) {
    const next = mode === "lie" ? "lie" : "strang";
    this.integrator = next;
    if (el.integratorVal) el.integratorVal.textContent = getIntegratorLabel(next);
    if (el.btnIntegrator) el.btnIntegrator.textContent = getIntegratorLabel(next);
  }

  applyPotentialPhase(scale) {
    const N = this.N;
    const NN = N * N;
    const factor = -this.dt * scale;
    for (let idx = 0; idx < NN; idx++) {
      const theta = factor * this.V[idx];
      const c = Math.cos(theta);
      const s = Math.sin(theta);
      const re = this.psiRe[idx];
      const im = this.psiIm[idx];
      this.psiRe[idx] = c * re - s * im;
      this.psiIm[idx] = s * re + c * im;
    }
  }

  applyKineticPhase(KRe, KIm) {
    const N = this.N;
    const NN = N * N;
    const kRe = KRe || this.KRe;
    const kIm = KIm || this.KIm;

    // Copy to WASM memory
    wasm.HEAPF64.set(this.psiRe, this.ptrRe >> 3);
    wasm.HEAPF64.set(this.psiIm, this.ptrIm >> 3);

    // Forward DST-I on rows (batch N, each length N)
    wasm._dst1_batch_forward(this.ptrRe, N, N);
    wasm._dst1_batch_forward(this.ptrIm, N, N);

    // Pull back, transpose to tmp
    this.psiRe.set(wasm.HEAPF64.subarray(this.ptrRe >> 3, (this.ptrRe >> 3) + NN));
    this.psiIm.set(wasm.HEAPF64.subarray(this.ptrIm >> 3, (this.ptrIm >> 3) + NN));

    transposeSquare(this.psiRe, this.tmpRe, N);
    transposeSquare(this.psiIm, this.tmpIm, N);

    // Push transposed to WASM
    wasm.HEAPF64.set(this.tmpRe, this.ptrRe >> 3);
    wasm.HEAPF64.set(this.tmpIm, this.ptrIm >> 3);

    // Forward DST-I on "y" (rows of transposed)
    wasm._dst1_batch_forward(this.ptrRe, N, N);
    wasm._dst1_batch_forward(this.ptrIm, N, N);

    // Pull back spectral arrays into tmpRe/tmpIm (still transposed layout, but that’s fine for phase multiply)
    this.tmpRe.set(wasm.HEAPF64.subarray(this.ptrRe >> 3, (this.ptrRe >> 3) + NN));
    this.tmpIm.set(wasm.HEAPF64.subarray(this.ptrIm >> 3, (this.ptrIm >> 3) + NN));

    // Multiply kinetic phase in spectral space.
    // Note: our tmp arrays correspond to (n,m) due to transpose.
    // K phase is symmetric in indices if L same, but we apply correctly by using swapped indexing:
    // tmp(n,m) *= exp(-i dt * 0.5*(kx_m^2 + ky_n^2)).
    for (let n = 0; n < N; n++) {
      for (let m = 0; m < N; m++) {
        const idx = m + N * n; // idx in transposed storage = (m,n) meaning original (n,m)
        const Kidx = n + N * m; // swap indices to match original (m,n)
        const c = kRe[Kidx];
        const s = kIm[Kidx];
        const re = this.tmpRe[idx];
        const im = this.tmpIm[idx];
        this.tmpRe[idx] = c * re - s * im;
        this.tmpIm[idx] = s * re + c * im;
      }
    }

    // Push back for inverse DSTs
    wasm.HEAPF64.set(this.tmpRe, this.ptrRe >> 3);
    wasm.HEAPF64.set(this.tmpIm, this.ptrIm >> 3);

    // Inverse DST-I on rows (y)
    wasm._dst1_batch_inverse(this.ptrRe, N, N);
    wasm._dst1_batch_inverse(this.ptrIm, N, N);

    // Pull back, transpose to psi
    this.tmpRe.set(wasm.HEAPF64.subarray(this.ptrRe >> 3, (this.ptrRe >> 3) + NN));
    this.tmpIm.set(wasm.HEAPF64.subarray(this.ptrIm >> 3, (this.ptrIm >> 3) + NN));

    transposeSquare(this.tmpRe, this.psiRe, N);
    transposeSquare(this.tmpIm, this.psiIm, N);

    // Push psi for inverse x
    wasm.HEAPF64.set(this.psiRe, this.ptrRe >> 3);
    wasm.HEAPF64.set(this.psiIm, this.ptrIm >> 3);

    // Inverse DST-I on rows (x)
    wasm._dst1_batch_inverse(this.ptrRe, N, N);
    wasm._dst1_batch_inverse(this.ptrIm, N, N);

    // Pull back to psi
    this.psiRe.set(wasm.HEAPF64.subarray(this.ptrRe >> 3, (this.ptrRe >> 3) + NN));
    this.psiIm.set(wasm.HEAPF64.subarray(this.ptrIm >> 3, (this.ptrIm >> 3) + NN));
  }

  finishStep() {
    const N = this.N;
    const NN = N * N;

    // Edge damping mask
    if (this.maskExtent > 0 && this.maskX && this.maskY) {
      for (let j = 0; j < N; j++) {
        const my = this.maskY[j];
        const row = j * N;
        for (let i = 0; i < N; i++) {
          const m = my * this.maskX[i];
          const idx = row + i;
          this.psiRe[idx] *= m;
          this.psiIm[idx] *= m;
        }
      }
    }

    this.t += this.dt;

    if (this.renormEnabled) {
      const norm = this.computeNorm();
      if (Number.isFinite(norm) && norm > 0) {
        const scale = 1.0 / Math.sqrt(norm);
        for (let idx = 0; idx < NN; idx++) {
          this.psiRe[idx] *= scale;
          this.psiIm[idx] *= scale;
        }
      }
    }
  }

  stepOnceStrang() {
    this.applyPotentialPhase(0.5);
    this.applyKineticPhase(this.KRe, this.KIm);
    this.applyPotentialPhase(0.5);
    this.finishStep();
  }

  stepOnceLie() {
    this.applyKineticPhase(this.KRe, this.KIm);
    this.applyPotentialPhase(1.0);
    this.finishStep();
  }

  // ---------- Splitting steps ----------
  stepOnce() {
    if (this.integrator === "lie") {
      this.stepOnceLie();
      return;
    }
    this.stepOnceStrang();
  }

  // ---------- Rendering ----------
  render(ctx) {
    const N = this.N;
    const w = el.canvas.width;
    const h = el.canvas.height;

    // Render at simulation resolution by drawing ImageData at N x N then scale.
    // For simplicity and speed, create an offscreen ImageData each frame.
    // If needed, cache ImageData and reuse.
    const img = ctx.createImageData(N, N);
    const data = img.data;

    // auto-scale magnitude to [0,1] using max
    let maxMag = 0;
    const range = this.getCurrentPotentialRange();
    const useDynamicRange = range.max <= range.min;
    let vMin = range.min;
    let vMax = range.max;
    const showMagnitude = this.visMode === "magnitude";
    if (useDynamicRange) {
      vMin = Infinity;
      vMax = -Infinity;
    }
    for (let idx = 0; idx < N * N; idx++) {
      const re = this.psiRe[idx];
      const im = this.psiIm[idx];
      const mag = Math.hypot(re, im);
      if (mag > maxMag) maxMag = mag;
      if (useDynamicRange) {
        const v = this.V[idx];
        if (v < vMin) vMin = v;
        if (v > vMax) vMax = v;
      }
    }
    const inv = maxMag > 0 ? 1 / maxMag : 1;
    const vRange = vMax - vMin;

    for (let j = 0; j < N; j++) {
      for (let i = 0; i < N; i++) {
        const idx = i + N * j;
        const re = this.psiRe[idx];
        const im = this.psiIm[idx];
        let r, g, b;
        const mag = Math.hypot(re, im);
        const magNorm = mag * inv;
        const v = Math.pow(clamp(mag, 0, 1), 0.6);
        if (showMagnitude) {
          r = v; g = v; b = v;
        } else {
          const phase = Math.atan2(im, re);
          const hue = (phase + Math.PI) / (2 * Math.PI); // [0,1)
          [r, g, b] = hsvToRgb(hue, 1.0, v);
        }

        // Potential overlay: transparent at Vmin, opaque white at Vmax.
        const vPot = this.V[idx];
        const a = vRange > 0 ? clamp((vPot - vMin) / vRange, 0, 1) : 0;
        if (a > 0) {
          r = r * (1 - a) + a;
          g = g * (1 - a) + a;
          b = b * (1 - a) + a;
        }

        const p = 4 * (i + N * (N - 1 - j)); // flip y for display
        data[p + 0] = Math.floor(255 * r);
        data[p + 1] = Math.floor(255 * g);
        data[p + 2] = Math.floor(255 * b);
        data[p + 3] = 255;
      }
    }

    // Draw simulation image scaled to canvas
    const off = document.createElement("canvas");
    off.width = N;
    off.height = N;
    const octx = off.getContext("2d", { willReadFrequently: false });
    octx.putImageData(img, 0, 0);

    ctx.imageSmoothingEnabled = false;
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, w, h);
    const size = Math.min(w, h);
    const offsetX = (w - size) * 0.5;
    const offsetY = (h - size) * 0.5;
    ctx.drawImage(off, 0, 0, N, N, offsetX, offsetY, size, size);

    // Overlay gaussian and velocity vector on recent slider changes.
    const now = performance.now();
    if (now < this.overlayUntil) {
      const L = this.L;
      const u = (this.x0 + L) / (2 * L);
      const v = (L - this.y0) / (2 * L);
      const cx = offsetX + u * size;
      const cy = offsetY + v * size;
      const r = Math.max(2, (this.sigma / (2 * L)) * size);

      const vMagMax = parseFloat(el.vMag?.max || "10");
      const speed = Math.max(0, this.vMag);
      const len = (vMagMax > 0 ? speed / vMagMax : 0) * (size * 0.25);
      const theta = (this.vAngleDeg * Math.PI) / 180;
      const dx = Math.cos(theta) * len;
      const dy = -Math.sin(theta) * len;

      ctx.save();
      ctx.strokeStyle = "rgba(255, 225, 0, 0.95)";
      ctx.fillStyle = "rgba(255, 225, 0, 0.95)";
      ctx.lineWidth = Math.max(1.5, size * 0.002);

      ctx.beginPath();
      ctx.arc(cx, cy, r, 0, Math.PI * 2);
      ctx.stroke();

      if (speed > 0) {
        ctx.beginPath();
        ctx.moveTo(cx, cy);
        ctx.lineTo(cx + dx, cy + dy);
        ctx.stroke();

        const ah = Math.max(6, size * 0.015);
        const angle = Math.atan2(dy, dx);
        ctx.beginPath();
        ctx.moveTo(cx + dx, cy + dy);
        ctx.lineTo(
          cx + dx - ah * Math.cos(angle - Math.PI / 6),
          cy + dy - ah * Math.sin(angle - Math.PI / 6)
        );
        ctx.lineTo(
          cx + dx - ah * Math.cos(angle + Math.PI / 6),
          cy + dy - ah * Math.sin(angle + Math.PI / 6)
        );
        ctx.closePath();
        ctx.fill();
      }
      ctx.restore();
    }

    // Hover brush outline
    if (this.hoverActive) {
      const L = this.L;
      const u = (this.hoverX + L) / (2 * L);
      const v = (L - this.hoverY) / (2 * L);
      const cx = offsetX + u * size;
      const cy = offsetY + v * size;
      const r = Math.max(2, ((this.brushSize * 0.5) / (2 * L)) * size);

      const { min: vMin, max: vMax } = this.getCurrentPotentialRange();
      const range = vMax - vMin;
      const t = range > 1e-6 ? clamp((this.brushValue - vMin) / range, 0, 1) : 1;
      const alpha = lerp(0.2, 1.0, t);

      ctx.save();
      ctx.strokeStyle = `rgba(255, 255, 255, ${alpha.toFixed(3)})`;
      ctx.lineWidth = Math.max(1.5, size * 0.002);
      ctx.beginPath();
      ctx.arc(cx, cy, r, 0, Math.PI * 2);
      ctx.stroke();
      ctx.restore();
    }
  }

  updateFps(now) {
    this.frameCount++;
    const dt = now - this.lastFpsTime;
    if (dt >= 500) {
      this.fps = (this.frameCount * 1000) / dt;
      this.frameCount = 0;
      this.lastFpsTime = now;
    }
  }

  computeNorm() {
    const N = this.N;
    let sum = 0;
    for (let idx = 0; idx < N * N; idx++) {
      const re = this.psiRe[idx];
      const im = this.psiIm[idx];
      sum += re * re + im * im;
    }
    return sum * this.dx * this.dx;
  }
}

function measureDstRoundTripGain(n) {
  const count = n;
  const bytes = count * 8;
  const ptr = wasm._malloc(bytes);
  const data = new Float64Array(count);
  for (let i = 0; i < count; i++) {
    data[i] = Math.sin((i + 1) * 0.123) + 0.3 * Math.cos((i + 1) * 0.57);
  }
  let inSum = 0;
  for (let i = 0; i < count; i++) inSum += data[i] * data[i];
  wasm.HEAPF64.set(data, ptr >> 3);
  wasm._dst1_batch_forward(ptr, n, 1);
  wasm._dst1_batch_inverse(ptr, n, 1);
  const out = wasm.HEAPF64.subarray(ptr >> 3, (ptr >> 3) + count);
  let outSum = 0;
  for (let i = 0; i < count; i++) outSum += out[i] * out[i];
  wasm._free(ptr);
  return inSum > 0 ? outSum / inSum : 0;
}

// Square transpose: out[j,i] = inp[i,j]
function transposeSquare(inp, out, N) {
  for (let j = 0; j < N; j++) {
    const rowOff = j * N;
    for (let i = 0; i < N; i++) {
      out[j + N * i] = inp[i + rowOff];
    }
  }
}

// ---------- Main ----------
const sim = new TDSE2D();
const ctx = el.canvas.getContext("2d", { alpha: false, desynchronized: true });

function getSquareViewport(canvas, rectOverride = null) {
  const rect = rectOverride || canvas.getBoundingClientRect();
  const size = Math.min(rect.width, rect.height);
  const offsetX = (rect.width - size) * 0.5;
  const offsetY = (rect.height - size) * 0.5;
  return { size, offsetX, offsetY, rect };
}

function applyUIConfig() {
  const rows = document.querySelectorAll("[data-ui]");
  rows.forEach((row) => {
    const key = row.getAttribute("data-ui");
    const show = Object.prototype.hasOwnProperty.call(UI_CONFIG, key)
      ? !!UI_CONFIG[key]
      : true;
    row.style.display = show ? "" : "none";
  });
}

function resetControlToDefault(ctrl) {
  if (!ctrl) return;
  const tag = ctrl.tagName;
  if (tag === "SELECT") {
    let idx = 0;
    for (let i = 0; i < ctrl.options.length; i++) {
      if (ctrl.options[i].defaultSelected) {
        idx = i;
        break;
      }
    }
    ctrl.selectedIndex = idx;
    return;
  }
  if (ctrl.type === "checkbox") {
    ctrl.checked = !!ctrl.defaultChecked;
    return;
  }
  if (ctrl.defaultValue !== undefined && ctrl.defaultValue !== "") {
    ctrl.value = ctrl.defaultValue;
  }
}

function attachUI() {
  function onChange() {
    sim.updateParamsFromUI();
  }
  for (const key of ["k", "L", "dt", "steps", "mask"]) {
    el[key].addEventListener("input", onChange);
  }
  el.renorm.addEventListener("change", onChange);

  POTENTIAL_PRESETS.forEach((preset, idx) => {
    const opt = document.createElement("option");
    opt.value = String(idx);
    opt.textContent = preset.name;
    el.preset.appendChild(opt);
  });
  if (DEFAULTS.presetIndex !== undefined) {
    const presetIdx = clamp(DEFAULTS.presetIndex, 0, POTENTIAL_PRESETS.length - 1);
    Array.from(el.preset.options).forEach((opt) => { opt.defaultSelected = false; });
    if (el.preset.options[presetIdx]) el.preset.options[presetIdx].defaultSelected = true;
    el.preset.value = String(presetIdx);
    sim.presetIndex = presetIdx;
    sim.preset = POTENTIAL_PRESETS[presetIdx];
  } else {
    el.preset.value = String(sim.presetIndex);
  }
  el.preset.addEventListener("change", async () => {
    await sim.setPreset(parseInt(el.preset.value, 10));
  });

  let savedSymbolicExamples = loadSavedSymbolicExamples();
  function getDefaultSymbolicExampleValue() {
    const idx = Number.isInteger(DEFAULTS.symbolicExampleIndex)
      ? clamp(DEFAULTS.symbolicExampleIndex, -1, SYMBOLIC_POTENTIAL_EXAMPLES.length - 1)
      : -1;
    return idx >= 0 ? `builtin:${idx}` : "custom";
  }
  function getSymbolicExpressionForOption(value) {
    if (typeof value !== "string") return null;
    if (value.startsWith("builtin:")) {
      const idx = parseInt(value.slice("builtin:".length), 10);
      if (Number.isInteger(idx) && idx >= 0 && idx < SYMBOLIC_POTENTIAL_EXAMPLES.length) {
        return SYMBOLIC_POTENTIAL_EXAMPLES[idx].expression;
      }
      return null;
    }
    if (value.startsWith("saved:")) {
      const idx = parseInt(value.slice("saved:".length), 10);
      if (Number.isInteger(idx) && idx >= 0 && idx < savedSymbolicExamples.length) {
        return savedSymbolicExamples[idx].expression;
      }
      return null;
    }
    return null;
  }
  function getSavedSymbolicIndexFromSelection() {
    const value = el.symbolicExprExample ? el.symbolicExprExample.value : "";
    if (typeof value !== "string" || !value.startsWith("saved:")) return -1;
    const idx = parseInt(value.slice("saved:".length), 10);
    if (!Number.isInteger(idx) || idx < 0 || idx >= savedSymbolicExamples.length) return -1;
    return idx;
  }
  function updateSymbolicActionButtons() {
    if (!el.btnDeleteSymbolic) return;
    el.btnDeleteSymbolic.disabled = getSavedSymbolicIndexFromSelection() < 0;
  }
  function syncExpressionFromDropdownSelection() {
    if (!el.symbolicExprExample || !el.symbolicExpr) return;
    const expression = getSymbolicExpressionForOption(el.symbolicExprExample.value);
    if (expression !== null) {
      el.symbolicExpr.value = expression;
    }
  }
  function rebuildSymbolicDropdown(selectedValue = null) {
    if (!el.symbolicExprExample) return;
    const defaultValue = getDefaultSymbolicExampleValue();
    el.symbolicExprExample.innerHTML = "";

    const customOption = document.createElement("option");
    customOption.value = "custom";
    customOption.textContent = "Custom expression";
    customOption.defaultSelected = defaultValue === "custom";
    el.symbolicExprExample.appendChild(customOption);

    SYMBOLIC_POTENTIAL_EXAMPLES.forEach((sample, idx) => {
      const opt = document.createElement("option");
      opt.value = `builtin:${idx}`;
      opt.textContent = sample.name;
      opt.defaultSelected = defaultValue === opt.value;
      el.symbolicExprExample.appendChild(opt);
    });

    savedSymbolicExamples.forEach((sample, idx) => {
      const opt = document.createElement("option");
      opt.value = `saved:${idx}`;
      opt.textContent = `${sample.name} (saved)`;
      opt.defaultSelected = false;
      el.symbolicExprExample.appendChild(opt);
    });

    const values = Array.from(el.symbolicExprExample.options).map((opt) => opt.value);
    const nextValue = selectedValue && values.includes(selectedValue)
      ? selectedValue
      : (values.includes(defaultValue) ? defaultValue : "custom");
    el.symbolicExprExample.value = nextValue;
    updateSymbolicActionButtons();
  }

  if (el.symbolicExprExample) {
    rebuildSymbolicDropdown();

    if (el.symbolicExpr && !(el.symbolicExpr.value || "").trim()) {
      const defaultExpression = getSymbolicExpressionForOption(el.symbolicExprExample.value);
      if (defaultExpression !== null) {
        el.symbolicExpr.value = defaultExpression;
        el.symbolicExpr.defaultValue = defaultExpression;
      }
    }

    el.symbolicExprExample.addEventListener("change", () => {
      syncExpressionFromDropdownSelection();
      updateSymbolicActionButtons();
    });
  }

  function evaluateSymbolicPotential() {
    const expression = el.symbolicExpr ? el.symbolicExpr.value : "";
    sim.applySymbolicExpression(expression);
  }
  if (el.btnEvalSymbolic) {
    el.btnEvalSymbolic.addEventListener("click", evaluateSymbolicPotential);
  }
  if (el.btnSaveSymbolic) {
    el.btnSaveSymbolic.addEventListener("click", () => {
      if (!el.symbolicExpr) return;
      const expression = el.symbolicExpr.value.trim();
      if (!expression) {
        el.presetInfo.textContent = "Cannot save an empty expression.";
        return;
      }

      const currentSelection = el.symbolicExprExample ? el.symbolicExprExample.value : "custom";
      let suggestedName = `Formula ${savedSymbolicExamples.length + 1}`;
      if (typeof currentSelection === "string" && currentSelection.startsWith("saved:")) {
        const idx = parseInt(currentSelection.slice("saved:".length), 10);
        if (Number.isInteger(idx) && idx >= 0 && idx < savedSymbolicExamples.length) {
          suggestedName = savedSymbolicExamples[idx].name;
        }
      }
      const typedName = prompt("Name for this saved expression:", suggestedName);
      if (typedName === null) return;
      const name = typedName.trim();
      if (!name) {
        el.presetInfo.textContent = "Save cancelled: name is empty.";
        return;
      }

      const existingIdx = savedSymbolicExamples.findIndex(
        (item) => item.name.toLowerCase() === name.toLowerCase()
      );
      if (existingIdx >= 0) {
        savedSymbolicExamples[existingIdx] = { name, expression };
        if (!saveSymbolicExamplesToStorage(savedSymbolicExamples)) {
          el.presetInfo.textContent = "Could not save to browser storage.";
          return;
        }
        rebuildSymbolicDropdown(`saved:${existingIdx}`);
        el.presetInfo.textContent = `Updated saved expression "${name}".`;
        return;
      }

      savedSymbolicExamples.push({ name, expression });
      if (!saveSymbolicExamplesToStorage(savedSymbolicExamples)) {
        savedSymbolicExamples.pop();
        el.presetInfo.textContent = "Could not save to browser storage.";
        return;
      }
      const newIdx = savedSymbolicExamples.length - 1;
      rebuildSymbolicDropdown(`saved:${newIdx}`);
      el.presetInfo.textContent = `Saved expression "${name}" to browser storage.`;
    });
  }
  if (el.btnDeleteSymbolic) {
    el.btnDeleteSymbolic.addEventListener("click", () => {
      const idx = getSavedSymbolicIndexFromSelection();
      if (idx < 0) {
        el.presetInfo.textContent = "Select a saved formula to delete.";
        return;
      }

      const item = savedSymbolicExamples[idx];
      const ok = confirm(`Delete saved formula "${item.name}"?`);
      if (!ok) return;

      const removed = savedSymbolicExamples.splice(idx, 1)[0];
      if (!saveSymbolicExamplesToStorage(savedSymbolicExamples)) {
        savedSymbolicExamples.splice(idx, 0, removed);
        el.presetInfo.textContent = "Could not delete from browser storage.";
        return;
      }

      rebuildSymbolicDropdown("custom");
      updateSymbolicActionButtons();
      el.presetInfo.textContent = `Deleted saved expression "${item.name}".`;
    });
    updateSymbolicActionButtons();
  }
  if (el.symbolicExpr) {
    el.symbolicExpr.addEventListener("keydown", (event) => {
      if ((event.ctrlKey || event.metaKey) && event.key === "Enter") {
        event.preventDefault();
        evaluateSymbolicPotential();
      }
    });
  }

  if (el.visMode) {
    el.visMode.value = sim.visMode;
    el.visMode.addEventListener("change", () => {
      sim.visMode = el.visMode.value;
    });
  }
  if (el.btnIntegrator) {
    sim.setIntegrator(sim.integrator);
    el.btnIntegrator.addEventListener("click", () => {
      const next = sim.integrator === "strang" ? "lie" : "strang";
      sim.setIntegrator(next);
      updateLabels(sim);
    });
  }

  el.btnToggle.addEventListener("click", () => {
    sim.running = !sim.running;
    el.btnToggle.textContent = sim.running ? "Pause" : "Start";
  });

  el.btnResetPsi.addEventListener("click", () => {
    sim.initPsiGaussian();
  });
  el.btnResetPotential.addEventListener("click", () => {
    sim.resetPotential();
  });
  el.btnResetAll.addEventListener("click", async () => {
    const controls = [
      el.k, el.L, el.dt, el.steps, el.mask,
      el.preset, el.visMode,
      el.symbolicExprExample, el.symbolicExpr,
      el.brushSize, el.brushHardness, el.brushValue,
      el.x0, el.y0, el.vMag, el.vAngle, el.sigma,
      el.renorm,
    ];
    controls.forEach(resetControlToDefault);
    syncExpressionFromDropdownSelection();

    await sim.updateParamsFromUI();
    if (el.preset) {
      await sim.setPreset(parseInt(el.preset.value, 10));
    }
    if (el.visMode) sim.visMode = el.visMode.value;
    sim.setIntegrator(DEFAULTS.integrator || "strang");
    sim.syncBrushValueRange();
    sim.updateBrushFromUI();
    sim.updateGaussianFromUI();
    sim.resetPotential();
  });

  for (const key of ["brushSize", "brushHardness", "brushValue"]) {
    el[key].addEventListener("input", () => sim.updateBrushFromUI());
  }
  for (const key of ["x0", "y0", "vMag", "vAngle", "sigma"]) {
    el[key].addEventListener("input", () => sim.updateGaussianFromUI());
  }

  // Painting interaction: map canvas position -> (x,y) in [-L,L]^2
  el.canvas.addEventListener("pointerdown", (e) => {
    el.canvas.setPointerCapture(e.pointerId);
    sim.painting = true;
    paintFromEvent(e);
  });
  el.canvas.addEventListener("pointerup", () => { sim.painting = false; });
  el.canvas.addEventListener("pointercancel", () => { sim.painting = false; });
  el.canvas.addEventListener("pointermove", (e) => {
    updateHoverFromEvent(e);
    if (sim.painting) paintFromEvent(e);
  });
  el.canvas.addEventListener("pointerenter", () => { sim.hoverActive = true; });
  el.canvas.addEventListener("pointerleave", () => { sim.hoverActive = false; });

  function paintFromEvent(e) {
    const { size, offsetX, offsetY, rect } = getSquareViewport(el.canvas);
    const px = e.clientX - rect.left - offsetX;
    const py = e.clientY - rect.top - offsetY;
    if (px < 0 || py < 0 || px > size || py > size) return;
    const u = px / size;   // [0,1]
    const v = py / size;   // [0,1]
    const x = -sim.L + (2 * sim.L) * u;
    const y =  sim.L - (2 * sim.L) * v; // invert
    sim.paintAt(x, y);
  }

  function updateHoverFromEvent(e) {
    const { size, offsetX, offsetY, rect } = getSquareViewport(el.canvas);
    const px = e.clientX - rect.left - offsetX;
    const py = e.clientY - rect.top - offsetY;
    if (px < 0 || py < 0 || px > size || py > size) return;
    const u = px / size;
    const v = py / size;
    sim.hoverX = -sim.L + (2 * sim.L) * u;
    sim.hoverY =  sim.L - (2 * sim.L) * v;
  }

  // hover updates handled in pointermove above

  window.addEventListener("resize", () => sim.resizeCanvas());
}

function updateStats() {
  const norm = sim.computeNorm();
  const preset = sim.preset ? sim.preset.name : "None";
  const vMin = sim.preset ? sim.preset.vMin : 0;
  const vMax = sim.preset ? sim.preset.vMax : 0;
  el.stats.innerHTML =
    //`t = ${sim.t.toFixed(3)}<br/>` +
    //`N = ${sim.N} (interior), dx = ${sim.dx.toFixed(4)}<br/>` +
    //`Δt = ${sim.dt.toFixed(4)}<br/>` +
    //`Preset = ${preset} (Vmin=${vMin.toFixed(2)}, Vmax=${vMax.toFixed(2)})<br/>` +
    `FPS = ${sim.fps.toFixed(1)}<br/>`;
    //`∥ψ∥² ≈ ${norm.toFixed(6)}<br/>` +
    //`DST round-trip gain ≈ ${sim.dstRoundTripGain.toFixed(6)}`;
}

function animate(now) {
  sim.updateFps(now);

  if (sim.running && !sim.loadingPreset) {
    const steps = Math.max(1, sim.stepsPerFrame | 0);
    for (let s = 0; s < steps; s++) sim.stepOnce();
  }

  sim.render(ctx);
  updateStats();

  requestAnimationFrame(animate);
}

// Boot
(async function main() {
  el.stats.textContent = "Loading WASM...";
  sim.resizeCanvas();
  applyUIConfig();
  attachUI();
  sim.syncRangesForL();
  sim.syncBrushValueRange();
  updateLabels(sim);

  await loadWasm();

  await sim.allocate();
  el.stats.textContent = "Ready.";
  requestAnimationFrame(animate);
})();
