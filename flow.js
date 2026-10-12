// Bubble-particle collision: steady axisymmetric flow past a rising bubble, particle motion and collision efficiency.
// Theoretical model only: no measured data and no data from any research centre are used.
// Pure functions, no DOM. Dimensionless variables: lengths in bubble radii, velocities in bubble rise speeds.
//
// Coordinates: spherical (r, theta) around the bubble centre, axisymmetric about the z axis.
// The liquid flows past the bubble in the +z direction (the bubble frame of a bubble rising in -z),
// so theta = pi is the front (upstream) stagnation point and theta = 0 is the rear.
// Stokes stream function psi: u_r = psi_theta / (r^2 sin theta), u_theta = -psi_r / (r sin theta),
// far field psi = r^2 sin^2(theta) / 2. zeta = omega r sin(theta) with E^2 psi = -zeta.

'use strict';

const GRAV = 9.81;

// ---------- drag correlations ----------
// Rigid sphere (bubble with an immobile surface, as in water with frother): Schiller and Naumann (1933), Re < 800.
const cdRigid = (Re) => 24 / Re * (1 + 0.15 * Math.pow(Re, 0.687));
// Clean spherical bubble (mobile surface): Mei, Klausner and Lawrence (1994).
const cdClean = (Re) => 16 / Re * (1 + 1 / (8 / Re + 0.5 * (1 + 3.315 / Math.sqrt(Re))));
// Separation angle behind a rigid sphere, measured from the front stagnation point, Clift, Grace and Weber (1978), 20 < Re < 400.
const sepAngleClift = (Re) => 180 - 42.5 * Math.pow(Math.log(Re / 20), 0.483);

// Terminal velocity from a force balance: |drho| g pi d^3/6 = Cd rhoL U^2 pi d^2 / 8.
function terminalVelocity(d, drho, rhoL, mu, cdFn) {
  const F = (U) => U * U - 4 * Math.abs(drho) * GRAV * d / (3 * cdFn(rhoL * U * d / mu) * rhoL);
  let lo = 1e-12, hi = 50;
  for (let i = 0; i < 200; i++) { const m = Math.sqrt(lo * hi); if (F(m) > 0) hi = m; else lo = m; }
  return Math.sqrt(lo * hi);
}

// ---------- collision efficiency correlations (interception only, relative to the bubble's projected area) ----------
const ecSutherland = (rp) => 3 * rp;                                     // potential flow, mobile surface
const ecGaudin = (rp) => 1.5 * rp * rp;                                  // Stokes flow, immobile surface
const ecYoonLuttrell = (rp, Re) => (1.5 + 4 * Math.pow(Re, 0.72) / 15) * rp * rp;   // intermediate Re, immobile surface
const ecWeberPaddock = (rp, Re) => 1.5 * (1 + (3 / 16) * Re / (1 + 0.249 * Math.pow(Re, 0.56))) * rp * rp;   // intermediate Re, immobile surface
// rp = particle radius / bubble radius

// ---------- analytic stream functions (dimensionless) ----------
const ANALYTIC = {
  stokesRigid: (r, th) => 0.5 * Math.sin(th) ** 2 * (r * r - 1.5 * r + 0.5 / r),
  hadamard: (r, th) => 0.5 * Math.sin(th) ** 2 * (r * r - r),
  potential: (r, th) => 0.5 * Math.sin(th) ** 2 * (r * r - 1 / r),
};

// Grid used by the page: 'normal' or 'fine'. A mobile (shear-free) surface has a thin boundary layer above
// Re of about 20 and needs a finer radial spacing.
function gridFor(Re, surface, quality) {
  const fine = quality === 'fine';
  if (surface === 'rigid') return fine ? { M: 120, dxi: 0.018 } : { M: 90, dxi: 0.025 };
  if (Re < 20) return fine ? { M: 120, dxi: 0.018 } : { M: 90, dxi: 0.025 };
  return fine ? { M: 120, dxi: 0.009 } : { M: 90, dxi: 0.0125 };
}

// ---------- Navier-Stokes solver ----------
// Stream function and vorticity on a grid uniform in xi = ln r and theta.
// surface: 'rigid' (no slip) or 'free' (shear-free, clean bubble).
function makeGrid(Re, opts) {
  opts = opts || {};
  const Rmax = opts.Rmax || (Re < 1 ? 300 : Re < 10 ? 120 : 60);
  const M = opts.M || 90;
  const N = opts.N || Math.ceil(Math.log(Rmax) / (opts.dxi || 0.025));
  return { N, M, dxi: Math.log(Rmax) / N, dth: Math.PI / M, Rmax };
}

function createState(Re, surface, opts) {
  const g = makeGrid(Re, opts);
  const { N, M, dxi, dth } = g;
  const W = M + 1;
  const st = {
    Re, surface, N, M, W, dxi, dth, Rmax: g.Rmax,
    psi: new Float64Array((N + 1) * W), zeta: new Float64Array((N + 1) * W),
    e2: new Float64Array(N + 1), em: new Float64Array(N + 1),
    s: new Float64Array(W), c: new Float64Array(W), cot: new Float64Array(W),
    iters: 0, change: Infinity, converged: false, scheme: (opts && opts.scheme) || 'deferred', lam: opts && opts.lam,
  };
  for (let i = 0; i <= N; i++) { st.e2[i] = Math.exp(2 * i * dxi); st.em[i] = Math.exp(-i * dxi); }
  for (let j = 0; j <= M; j++) {
    const t = j * dth;
    st.s[j] = Math.sin(t); st.c[j] = Math.cos(t);
    st.cot[j] = j > 0 && j < M ? st.c[j] / st.s[j] : 0;
  }
  // start from the Stokes solution of the right kind
  for (let i = 0; i <= N; i++) {
    const r = Math.exp(i * dxi);
    for (let j = 0; j <= M; j++) {
      const s2 = st.s[j] * st.s[j];
      const k = i * W + j;
      if (surface === 'rigid') { st.psi[k] = 0.5 * s2 * (r * r - 1.5 * r + 0.5 / r); st.zeta[k] = -1.5 * s2 / r; }
      else { st.psi[k] = 0.5 * s2 * (r * r - r); st.zeta[k] = -s2 / r; }
    }
  }
  for (let j = 0; j <= M; j++) { st.psi[N * W + j] = 0.5 * st.e2[N] * st.s[j] * st.s[j]; st.zeta[N * W + j] = 0; }
  return st;
}

// One block of iterations. Returns the largest vorticity change in the last sweep, relative to max |zeta|.
function iterate(st, nIter) {
  const { Re, N, M, W, dxi, dth, psi, zeta, e2, em, s, cot } = st;
  const ix2 = 1 / (dxi * dxi), it2 = 1 / (dth * dth), ix = 1 / (2 * dxi), itd = 1 / (2 * dth);
  const wPsi = 1.8;
  const lam = st.lam || 0.8;
  const lamW = lam;
  let rel = Infinity;
  for (let n = 0; n < nIter; n++) {
    // stream function: one SOR sweep
    for (let i = 1; i < N; i++) {
      const aE = ix2 - ix, aWc = ix2 + ix, aP = 2 * ix2 + 2 * it2;
      for (let j = 1; j < M; j++) {
        const k = i * W + j;
        const aN = it2 - cot[j] * itd, aS = it2 + cot[j] * itd;
        const v = (aE * psi[k + W] + aWc * psi[k - W] + aN * psi[k + 1] + aS * psi[k - 1] + e2[i] * zeta[k]) / aP;
        psi[k] += wPsi * (v - psi[k]);
      }
    }
    // vorticity on the bubble surface
    for (let j = 1; j < M; j++) {
      let zw;
      if (st.surface === 'rigid') zw = -(8 * psi[W + j] - psi[2 * W + j]) * ix2 / 2;      // no slip: zeta = -psi_xixi
      else zw = -2 * (4 * psi[W + j] - psi[2 * W + j]) * ix;                              // shear free: zeta = -2 psi_xi
      zeta[j] += lamW * (zw - zeta[j]);
    }
    // outer boundary: no vorticity upstream, zero gradient in the wake
    for (let j = 1; j < M; j++) zeta[N * W + j] = j < M / 2 ? zeta[(N - 1) * W + j] : 0;
    // vorticity transport: one relaxation sweep.
    // scheme 'hybrid': central differences, first-order upwind where the cell Peclet number exceeds 2.
    // scheme 'deferred': upwind coefficients in the matrix plus the difference (central - upwind) as a
    // lagged source, so the converged solution satisfies the second-order central scheme everywhere.
    const deferred = st.scheme === 'deferred';
    let dmax = 0, zmax = 0;
    for (let i = 1; i < N; i++) {
      for (let j = 1; j < M; j++) {
        const k = i * W + j;
        const kk = 0.5 * Re * em[i] / s[j];
        const pT = (psi[k + 1] - psi[k - 1]) * itd, pX = (psi[k + W] - psi[k - W]) * ix;
        const A = -1 - kk * pT, B = -cot[j] + kk * pX, C = 2 * kk * (pT - cot[j] * pX);
        const zP = zeta[k], zE = zeta[k + W], zWv = zeta[k - W], zN = zeta[k + 1], zS = zeta[k - 1];
        let cE, cWw, cN, cS, diag, src = 0;
        const upX = deferred || Math.abs(A) * dxi > 2, upT = deferred || Math.abs(B) * dth > 2;
        if (!upX) { cE = ix2 + A * ix; cWw = ix2 - A * ix; diag = 2 * ix2; }
        else {
          cE = ix2 + (A > 0 ? A / dxi : 0); cWw = ix2 + (A < 0 ? -A / dxi : 0); diag = 2 * ix2 + Math.abs(A) / dxi;
          if (deferred) {
            const central = A * (zE - zWv) * ix;
            const upwind = A > 0 ? A * (zE - zP) / dxi : A * (zP - zWv) / dxi;
            src += central - upwind;
          }
        }
        if (!upT) { cN = it2 + B * itd; cS = it2 - B * itd; diag += 2 * it2; }
        else {
          cN = it2 + (B > 0 ? B / dth : 0); cS = it2 + (B < 0 ? -B / dth : 0); diag += 2 * it2 + Math.abs(B) / dth;
          if (deferred) {
            const central = B * (zN - zS) * itd;
            const upwind = B > 0 ? B * (zN - zP) / dth : B * (zP - zS) / dth;
            src += central - upwind;
          }
        }
        let rhs = cE * zE + cWw * zWv + cN * zN + cS * zS + src;
        if (C < 0) diag -= C; else rhs += C * zP;
        const nz = zP + lam * (rhs / diag - zP);
        const d = Math.abs(nz - zP);
        if (d > dmax) dmax = d;
        zeta[k] = nz;
        if (Math.abs(nz) > zmax) zmax = Math.abs(nz);
      }
    }
    st.iters++;
    rel = zmax > 0 ? dmax / zmax : dmax;
    if (!isFinite(zmax) || !isFinite(dmax) || zmax > 1e8) { st.diverged = true; st.change = NaN; return NaN; }
  }
  st.change = rel;
  return rel;
}

// Default vorticity relaxation 0.8; if the iteration diverges it restarts with half the relaxation.
function solveFlow(Re, surface, opts) {
  opts = Object.assign({}, opts || {});
  let lam = opts.lam || 0.8;
  for (;;) {
    const st = createState(Re, surface, Object.assign({}, opts, { lam }));
    const tol = opts.tol || 1e-7, maxIt = opts.maxIt || 60000;
    while (st.iters < maxIt) {
      const rel = iterate(st, 200);
      if (st.diverged) break;
      if (rel < tol) { st.converged = true; break; }
    }
    if (!st.diverged || lam < 0.1) return st;
    lam *= 0.5;
  }
}

// ---------- quantities from a solution ----------
function wallData(st) {
  const { N, M, W, dxi, dth, psi, zeta, s, c } = st;
  const zw = new Float64Array(M + 1), zx = new Float64Array(M + 1), ut = new Float64Array(M + 1);
  for (let j = 0; j <= M; j++) {
    zw[j] = zeta[j];
    zx[j] = (-3 * zeta[j] + 4 * zeta[W + j] - zeta[2 * W + j]) / (2 * dxi);
    ut[j] = j > 0 && j < M ? -((4 * psi[W + j] - psi[2 * W + j]) / (2 * dxi)) / s[j] : 0;   // u_theta at the wall
  }
  return { zw, zx, ut };
}

function dragCoefficient(st) {
  const { Re, M, dth, s, c } = st;
  const { zw, zx, ut } = wallData(st);
  const trap = (f) => { let a = 0; for (let j = 0; j <= M; j++) a += (j === 0 || j === M ? 0.5 : 1) * f(j); return a * dth; };
  if (st.surface === 'rigid') {
    const pressure = (8 / Re) * 0.5 * trap((j) => zx[j] * s[j]);
    const friction = -(8 / Re) * trap((j) => zw[j] * s[j]);
    return { cd: pressure + friction, pressure, friction };
  }
  // shear-free surface: pressure from the tangential momentum balance (with convection), plus normal viscous stress
  const dut = (j) => (j === 0 ? (ut[1] - ut[0]) / dth : j === M ? (ut[M] - ut[M - 1]) / dth : (ut[j + 1] - ut[j - 1]) / (2 * dth));
  const uts = (j) => ut[j] * s[j];
  const duts = (j) => (j === 0 ? (uts(1) - uts(0)) / dth : j === M ? (uts(M) - uts(M - 1)) / dth : (uts(j + 1) - uts(j - 1)) / (2 * dth));
  const conv = 2 * trap((j) => -ut[j] * dut(j) * s[j] * s[j]);
  const visc = (4 / Re) * trap((j) => zx[j] * s[j]);
  const normal = -(16 / Re) * trap((j) => duts(j) * c[j]);
  return { cd: conv + visc + normal, pressure: conv + visc, friction: normal };
}

// Separation angle from the front stagnation point in degrees (rigid surface), or null if the flow stays attached.
function separationAngle(st) {
  if (st.surface !== 'rigid') return null;
  const { M, dth } = st;
  const { zw } = wallData(st);
  // the front part of the wall (theta near pi) has zeta < 0; a recirculation region at the rear has zeta > 0
  for (let j = 1; j < M - 1; j++) {
    if (zw[j] > 0 && zw[j + 1] <= 0) {
      const t = (j + zw[j] / (zw[j] - zw[j + 1])) * dth;
      return 180 - t * 180 / Math.PI;
    }
  }
  return null;
}

// Length of the recirculation region behind the bubble, measured from the rear surface, in bubble diameters.
function wakeLength(st) {
  const { N, W, dxi, psi } = st;
  // psi < 0 just off the rear axis (j = 1) marks reversed flow
  let last = -1;
  for (let i = 1; i < N; i++) if (psi[i * W + 1] < 0) last = i; else if (last >= 0) break;
  if (last < 0) return 0;
  const r = Math.exp((last + 0.5) * dxi);
  return (r - 1) / 2;
}

// ---------- evaluating the field anywhere ----------
// Returns a function (r, theta) -> psi using cubic interpolation in xi and theta; a near-wall Taylor
// expansion is used inside the first cell, where psi varies as (r - 1)^2 (no slip) or (r - 1) (free).
function psiEvaluator(st) {
  if (st.analytic) return ANALYTIC[st.analytic];
  const { N, M, W, dxi, dth, psi } = st;
  const { zw, zx } = wallData(st);
  const cub = (p0, p1, p2, p3, t) => {           // Catmull-Rom style cubic through 4 equally spaced points, t in [0,1] between p1 and p2
    return p1 + 0.5 * t * (p2 - p0 + t * (2 * p0 - 5 * p1 + 4 * p2 - p3 + t * (3 * (p1 - p2) + p3 - p0)));
  };
  const atI = (i, th) => {                       // psi at grid row i, interpolated in theta
    const x = th / dth; let j = Math.floor(x); if (j < 0) j = 0; if (j > M - 1) j = M - 1;
    const t = x - j;
    const g = (jj) => (jj < 0 ? psi[i * W - jj] : jj > M ? psi[i * W + 2 * M - jj] : psi[i * W + jj]);   // psi ~ sin^2(theta) is even about the axis
    return cub(g(j - 1), g(j), g(j + 1), g(j + 2), t);
  };
  const wallTh = (arr, th) => { const x = th / dth; let j = Math.floor(x); if (j < 0) j = 0; if (j > M - 1) j = M - 1; const t = x - j; return arr[j] * (1 - t) + arr[j + 1] * t; };
  return (r, th) => {
    if (r <= 1) return 0;
    if (r >= st.Rmax) return 0.5 * r * r * Math.sin(th) ** 2;
    const x = Math.log(r) / dxi;
    if (x < 1) {
      // Taylor expansion about the wall in xi: psi = a1 xi + a2 xi^2 + a3 xi^3, matched to the wall data and to row 1
      const h = Math.log(r);
      if (st.surface === 'rigid') {
        // psi_xixi(0) = -zeta_w ; psi_xixixi(0) from zeta_xi: psi = -zeta_w h^2/2 + b h^3 with b fitted to row 1
        const a2 = -wallTh(zw, th) / 2;
        const p1 = atI(1, th);
        const b = (p1 - a2 * dxi * dxi) / (dxi * dxi * dxi);
        return a2 * h * h + b * h * h * h;
      }
      // free surface: psi = a1 h + a2 h^2 with a1 from rows 1 and 2 (one-sided), a2 fitted to row 1
      const p1 = atI(1, th), p2 = atI(2, th);
      const a1 = (4 * p1 - p2) / (2 * dxi);
      const a2 = (p1 - a1 * dxi) / (dxi * dxi);
      return a1 * h + a2 * h * h;
    }
    let i = Math.floor(x); if (i > N - 2) i = N - 2;
    const t = x - i;
    const row = (ii) => atI(Math.min(Math.max(ii, 0), N), th);
    if (i < 1) return row(1) * t;   // not reached
    return cub(row(i - 1), row(i), row(i + 1), row(i + 2), t);
  };
}

// ---------- velocity field on the grid ----------
// Velocity (u_x, u_z) and convective acceleration (a_x, a_z) = (u . grad) u at every grid node, from psi by
// second-order differences. x is the distance from the axis, z points along the flow (downstream).
function nodeFields(st) {
  const { N, M, W, dxi, dth, psi, s, c } = st;
  const n = (N + 1) * W;
  const ux = new Float64Array(n), uz = new Float64Array(n), ax = new Float64Array(n), az = new Float64Array(n);
  const P = (i, j) => (j < 0 ? psi[i * W - j] : j > M ? psi[i * W + 2 * M - j] : psi[i * W + j]);   // even about the axis
  for (let i = 0; i <= N; i++) {
    const r = Math.exp(i * dxi), r2 = r * r;
    for (let j = 0; j <= M; j++) {
      const k = i * W + j;
      let ur, ut;
      if (j === 0 || j === M) {
        // on the axis psi ~ c theta^2: u_r = psi_thth / r^2 (sign by side), u_theta = 0
        const ptt = (P(i, j + 1) - 2 * P(i, j) + P(i, j - 1)) / (dth * dth);
        ur = (j === 0 ? 1 : -1) * ptt / r2; ut = 0;
      } else {
        const pT = (P(i, j + 1) - P(i, j - 1)) / (2 * dth);
        let pX;
        if (i === 0) pX = (-3 * psi[k] + 4 * psi[k + W] - psi[k + 2 * W]) / (2 * dxi);
        else if (i === N) pX = (3 * psi[k] - 4 * psi[k - W] + psi[k - 2 * W]) / (2 * dxi);
        else pX = (psi[k + W] - psi[k - W]) / (2 * dxi);
        ur = pT / (r2 * s[j]); ut = -pX / (r2 * s[j]);
      }
      if (i === 0) { ur = 0; if (st.surface === 'rigid') ut = 0; }
      ux[k] = ur * s[j] + ut * c[j];
      uz[k] = ur * c[j] - ut * s[j];
    }
  }
  // convective acceleration from the Cartesian components: d/dx = sin d/dr + (cos/r) d/dth, d/dz = cos d/dr - (sin/r) d/dth
  const D = (f, i, j) => {
    const k = i * W + j;
    const g = (jj) => (jj < 0 ? f[i * W - jj] * 1 : jj > M ? f[i * W + 2 * M - jj] : f[i * W + jj]);
    let fX;
    if (i === 0) fX = (-3 * f[k] + 4 * f[k + W] - f[k + 2 * W]) / (2 * dxi);
    else if (i === N) fX = (3 * f[k] - 4 * f[k - W] + f[k - 2 * W]) / (2 * dxi);
    else fX = (f[k + W] - f[k - W]) / (2 * dxi);
    let fT;
    if (j === 0 || j === M) fT = 0;
    else fT = (g(j + 1) - g(j - 1)) / (2 * dth);
    const r = Math.exp(i * dxi);
    const fr = fX / r;
    return [s[j] * fr + c[j] * fT / r, c[j] * fr - s[j] * fT / r];
  };
  for (let i = 0; i <= N; i++) {
    for (let j = 0; j <= M; j++) {
      const k = i * W + j;
      const gx = D(ux, i, j), gz = D(uz, i, j);
      // u_x is odd across the axis: its theta derivative at the axis is not zero, use one-sided values there
      if (j === 0 || j === M) {
        const jj = j === 0 ? 1 : M - 1;
        const g1 = D(ux, i, jj), g2 = D(uz, i, jj);
        gx[0] = g1[0]; gx[1] = g1[1]; gz[0] = g2[0]; gz[1] = g2[1];
      }
      ax[k] = ux[k] * gx[0] + uz[k] * gx[1];
      az[k] = ux[k] * gz[0] + uz[k] * gz[1];
    }
  }
  return { ux, uz, ax, az };
}

// Returns f(x, z) -> [u_x, u_z, a_x, a_z]. The velocity is taken from derivatives of the stream-function
// interpolant, so it is exactly the velocity whose streamlines are the contours of psiEvaluator: a tracer then
// keeps its value of psi, and trajectories agree with the stream-function results (bilinear interpolation of
// nodal velocities drifts across streamlines by up to about 1 % over a long path). Near the surface the
// interpolant has the correct form (psi ~ (r - 1)^2 on a no-slip surface). The acceleration, which only enters
// through the small pressure-gradient and added-mass term, is interpolated bilinearly from the nodes.
// The sign of u_x and a_x follows the sign of x, so the field can be sampled on both sides of the axis.
function flowSampler(st) {
  const F = nodeFields(st);
  const ev = psiEvaluator(st);
  const { N, M, W, dxi, dth, Rmax } = st;
  const thMin = 2e-3;
  return (x, z) => {
    const sx = x < 0 ? -1 : 1, X = Math.abs(x);
    const r = Math.hypot(X, z);
    if (r >= Rmax) return [0, 1, 0, 0];
    if (r <= 1) return [0, 0, 0, 0];
    const th = Math.atan2(X, z);
    let fi = Math.log(r) / dxi, fj = th / dth;
    let i = Math.floor(fi), j = Math.floor(fj);
    if (i > N - 1) i = N - 1; if (j > M - 1) j = M - 1; if (j < 0) j = 0;
    const ti = fi - i, tj = fj - j;
    const k00 = i * W + j, k10 = k00 + W, k01 = k00 + 1, k11 = k10 + 1;
    const bl = (f) => (f[k00] * (1 - ti) + f[k10] * ti) * (1 - tj) + (f[k01] * (1 - ti) + f[k11] * ti) * tj;
    // very close to the axis sin(theta) -> 0: evaluate at thMin; u_z is even about the axis and u_x is odd and
    // linear in theta there, so u_x is scaled back to the true angle
    const t = Math.min(Math.max(th, thMin), Math.PI - thMin);
    const sn = Math.sin(t), cs = Math.cos(t);
    const h = Math.min(1e-4 * r, 0.25 * (r - 1)), d = 1e-4;
    const pr = (ev(r + h, t) - ev(r - h, t)) / (2 * h);
    const pt = (ev(r, t + d) - ev(r, t - d)) / (2 * d);
    const ur = pt / (r * r * sn), ut = -pr / (r * sn);
    let ux = ur * sn + ut * cs;
    const uz = ur * cs - ut * sn;
    if (t !== th) ux *= Math.sin(th) / sn;
    return [sx * ux, uz, sx * bl(F.ax), bl(F.az)];
  };
}

// ---------- particles ----------
// Particle properties in the dimensionless frame of the bubble (lengths in bubble radii a, velocities in the
// rise speed U, time in a / U). Equation of motion (Maxey and Riley without the history and lift forces,
// added-mass coefficient 1/2):
//   dv/dt = beta Du/Dt + (u - v + ws e_z) / St
//   beta = 1.5 / (rho_p / rho_f + 0.5),  St = (rho_p + 0.5 rho_f) d_p^2 U / (18 mu f a),  f = 1 + 0.15 Re_p^0.687
// ws is the terminal settling speed of the particle in still liquid divided by U (e_z points down, along the
// liquid flow past the bubble).
function particleParams(bub, par) {
  const { db, U, rhoL, mu } = bub;                       // SI
  const { dp, rhop } = par;
  const drho = rhop - rhoL;
  const wsDim = drho === 0 ? 0 : Math.sign(drho) * terminalVelocity(dp, Math.abs(drho), rhoL, mu, cdRigid);
  const Rep = rhoL * Math.abs(wsDim) * dp / mu;
  const f = 1 + 0.15 * Math.pow(Rep, 0.687);
  const a = db / 2;
  return {
    rp: dp / db,                                           // particle radius / bubble radius
    ws: wsDim / U,
    wsDim, Rep,
    St: (rhop + 0.5 * rhoL) * dp * dp * U / (18 * mu * f * a),
    beta: 1.5 / (rhop / rhoL + 0.5),
  };
}

// Collision efficiency of inertia-free particles (interception and settling), exact for the computed flow:
// the particle velocity u + ws e_z has the stream function Psi_p = psi + ws r^2 sin^2(theta) / 2, so the
// particles that reach the collision sphere r = 1 + rp are those inside the particle streamline that grazes
// it. E is relative to the particle flux through the bubble's projected area far upstream.
function collisionStream(evalPsi, rp, ws) {
  const r = 1 + rp;
  let best = 0, prev = -Infinity, thBest = Math.PI;
  const n = 3600;
  for (let k = 1; k < n; k++) {
    const th = Math.PI * (1 - k / n);                     // from the front stagnation point backwards
    const v = evalPsi(r, th) + 0.5 * ws * r * r * Math.sin(th) ** 2;
    if (v > best) { best = v; thBest = th; }
    if (v < prev && v < best) break;                       // past the first maximum: the grazing point
    prev = v;
  }
  return { E: 2 * best / (1 + ws), thetaGraze: thBest };
}

// One particle trajectory with inertia. Starts upstream at (x0, z0) with the local liquid velocity plus the
// settling speed. Returns { hit, path }; with opts.noPath the path is not stored (for the bisection), otherwise
// points about 0.004 radii apart.
// For St >= 0.02 the equation is integrated with a fourth-order Runge-Kutta method. Below that it is stiff (the
// particle relaxes to the liquid velocity much faster than the liquid velocity changes), and an exponential
// integrator is used instead: the drag term is integrated exactly, with the rest of the right-hand side taken
// as linear in time over the step (second order), and the step follows the distance the particle moves.
// opts.method = 'rk4' or 'exp' forces one of the two.
function trajectory(sample, evalPsi, x0, pp, opts) {
  opts = opts || {};
  const z0 = opts.z0 || -8, zEnd = opts.zEnd || 4, xMax = opts.xMax || 8, keep = !opts.noPath;
  const { rp, ws, St, beta } = pp;
  const method = opts.method || (St < 0.02 ? 'exp' : 'rk4');
  let x = x0, z = z0;
  const u0 = sample(x, z);
  let vx = u0[0], vz = u0[1] + ws;
  const path = [[x, z]];
  let lx = x, lz = z;
  const maxSteps = opts.maxSteps || 400000;
  const done = (hit) => { if (keep) path.push([x, z]); return { hit, path }; };
  if (method === 'rk4') {
    const rhs = (X, Z, VX, VZ) => {
      const q = sample(X, Z);
      return [VX, VZ, beta * q[2] + (q[0] - VX) / St, beta * q[3] + (q[1] - VZ + ws) / St];
    };
    for (let n = 0; n < maxSteps; n++) {
      const r = Math.hypot(x, z);
      if (r <= 1 + rp) return done(true);
      if (z > zEnd || x > xMax) return done(false);
      const h = Math.min(0.1, Math.max(2e-4, 0.2 * (r - 1 - rp)), 0.25 * St);
      const k1 = rhs(x, z, vx, vz);
      const k2 = rhs(x + 0.5 * h * k1[0], z + 0.5 * h * k1[1], vx + 0.5 * h * k1[2], vz + 0.5 * h * k1[3]);
      const k3 = rhs(x + 0.5 * h * k2[0], z + 0.5 * h * k2[1], vx + 0.5 * h * k2[2], vz + 0.5 * h * k2[3]);
      const k4 = rhs(x + h * k3[0], z + h * k3[1], vx + h * k3[2], vz + h * k3[3]);
      x += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
      z += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
      vx += h / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);
      vz += h / 6 * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3]);
      if (x < 0) { x = -x; vx = -vx; }
      if (keep && (x - lx) ** 2 + (z - lz) ** 2 > 1.6e-5) { path.push([x, z]); lx = x; lz = z; }
    }
    return { hit: false, path };
  }
  // exponential integrator: dx/dt = v, dv/dt = (g(x) - v) / St with g = u + ws e_z + St beta Du/Dt
  const g = (X, Z) => { const q = sample(X, Z); return [q[0] + St * beta * q[2], q[1] + ws + St * beta * q[3]]; };
  let g0 = g(x, z);
  for (let n = 0; n < maxSteps; n++) {
    const r = Math.hypot(x, z);
    if (r <= 1 + rp) return done(true);
    if (z > zEnd || x > xMax) return done(false);
    const speed = Math.max(Math.hypot(vx, vz), 1e-6);
    const h = Math.min(0.05, 0.05 * Math.max(r - 1 - rp, 2e-4) / speed);
    const E = Math.exp(-h / St), S = St * (1 - E);              // S = integral of exp(-t/St) over the step
    // predictor with g frozen, then g taken as linear in time between the two ends of the step
    const g1 = g(x + g0[0] * h + (vx - g0[0]) * S, z + g0[1] * h + (vz - g0[1]) * S);
    const sx = (g1[0] - g0[0]) / h, sz = (g1[1] - g0[1]) / h;
    // exact solution for g = g0 + s t:  v = g0 + s (t - St) + (v0 - g0 + St s) exp(-t/St)
    const cx = vx - g0[0] + St * sx, cz = vz - g0[1] + St * sz;
    x += g0[0] * h + sx * h * (0.5 * h - St) + cx * S;
    z += g0[1] * h + sz * h * (0.5 * h - St) + cz * S;
    vx = g0[0] + sx * (h - St) + cx * E;
    vz = g0[1] + sz * (h - St) + cz * E;
    if (x < 0) { x = -x; vx = -vx; }
    g0 = g(x, z);
    if (keep && (x - lx) ** 2 + (z - lz) ** 2 > 1.6e-5) { path.push([x, z]); lx = x; lz = z; }
  }
  return { hit: false, path };
}

// Inertia-free trajectory (follows u + ws e_z), for drawing and for checking the integration.
function tracer(sample, x0, ws, rp, opts) {
  opts = opts || {};
  const z0 = opts.z0 || -8, zEnd = opts.zEnd || 4, keep = !opts.noPath;
  let x = x0, z = z0;
  const path = [[x, z]];
  let lx = x, lz = z;
  const v = (X, Z) => { const q = sample(X, Z); return [q[0], q[1] + ws]; };
  for (let n = 0; n < 400000; n++) {
    const r = Math.hypot(x, z);
    if (r <= 1 + rp) { if (keep) path.push([x, z]); return { hit: true, path }; }
    if (z > zEnd || x > 8) { if (keep) path.push([x, z]); return { hit: false, path }; }
    const h = Math.min(0.05, Math.max(2e-4, 0.2 * (r - 1 - rp)));
    const k1 = v(x, z), k2 = v(x + 0.5 * h * k1[0], z + 0.5 * h * k1[1]);
    const k3 = v(x + 0.5 * h * k2[0], z + 0.5 * h * k2[1]), k4 = v(x + h * k3[0], z + h * k3[1]);
    x += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    z += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
    if (x < 0) x = -x;
    if (keep && (x - lx) ** 2 + (z - lz) ** 2 > 1.6e-5) { path.push([x, z]); lx = x; lz = z; }
  }
  return { hit: false, path };
}

// Starting offset (distance from the axis at z = z0) of the particle streamline Psi_p = E (1 + ws) / 2, i.e. of
// the grazing particle without inertia when E is its collision efficiency.
function streamOffset(evalPsi, ws, E, z0) {
  const Psi = (x) => evalPsi(Math.hypot(x, z0), Math.atan2(x, z0)) + 0.5 * ws * x * x;
  const target = 0.5 * E * (1 + ws);
  let lo = 0, hi = 3;
  for (let k = 0; k < 60; k++) { const m = 0.5 * (lo + hi); if (Psi(m) < target) lo = m; else hi = m; }
  return 0.5 * (lo + hi);
}

// Collision efficiency with inertia: bisection on the starting offset x0 at z = z0 upstream. Upstream of z0 the
// flow is nearly uniform, so the particle moves with the liquid plus its settling speed, and E follows from the
// particle stream function at the start. The search is bracketed around the grazing offset without inertia.
// Written as a generator that yields after each trajectory, so a page can keep responding; collisionInertial
// runs it to the end.
function* collisionInertialSteps(sample, evalPsi, pp, opts) {
  opts = opts || {};
  const z0 = opts.z0 || -8, rel = opts.rel || 2e-4;
  const Psi = (x) => evalPsi(Math.hypot(x, z0), Math.atan2(x, z0)) + 0.5 * pp.ws * x * x;
  let n = 0;
  const hits = (x) => { n++; return trajectory(sample, evalPsi, x, pp, { z0, noPath: true, method: opts.method }).hit; };
  const result = (x0, extra) => Object.assign({ E: 2 * Psi(x0) / (1 + pp.ws), x0, trajectories: n }, extra || {});
  const guess = streamOffset(evalPsi, pp.ws, collisionStream(evalPsi, pp.rp, pp.ws).E, z0);
  let lo = 0.8 * guess, hi = 1.3 * guess, hiTested = false;
  // widen the bracket until lo hits and hi misses
  for (;;) { const h = hits(lo); yield n; if (h) break; hi = lo; hiTested = true; lo *= 0.6; if (lo < 1e-9) return result(0); }
  if (!hiTested) for (;;) { const h = hits(hi); yield n; if (!h) break; lo = hi; hi *= 1.4; if (hi > 3) return result(3, { saturated: true }); }
  while (hi - lo > rel * hi) {
    const m = 0.5 * (lo + hi);
    if (hits(m)) lo = m; else hi = m;
    yield n;
  }
  return result(0.5 * (lo + hi), { lo, hi });     // lo hits, hi misses
}
function collisionInertial(sample, evalPsi, pp, opts) {
  const it = collisionInertialSteps(sample, evalPsi, pp, opts);
  for (;;) { const s = it.next(); if (s.done) return s.value; }
}

if (typeof module !== 'undefined') {
  module.exports = {
    GRAV, cdRigid, cdClean, sepAngleClift, terminalVelocity, ecSutherland, ecGaudin, ecYoonLuttrell, ecWeberPaddock, ANALYTIC,
    gridFor, makeGrid, createState, iterate, solveFlow, wallData, dragCoefficient, separationAngle, wakeLength,
    psiEvaluator, nodeFields, flowSampler, particleParams, collisionStream, streamOffset, trajectory, tracer,
    collisionInertialSteps, collisionInertial,
  };
}
