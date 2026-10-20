// Live two-dimensional flow past a bubble held in place, with an optional grid of rods upstream:
// lattice Boltzmann method, D2Q9 lattice, two-relaxation-time collision (Ginzburg, Verhaeghe and d'Humieres 2008)
// with the "magic" parameter 1/4, halfway bounce-back at solid surfaces (no slip), a velocity inlet after Zou and He (1997), and an
// outlet of equilibrium populations at density 1 with the velocity of the cell upstream, behind a sponge layer. Inertial particles follow the Maxey-Riley equation without history and lift forces.
// Optional Smagorinsky subgrid model for large-eddy simulation (Yu, Girimaji and Luo 2005).
// Turbulence statistics: mean velocity, turbulent kinetic energy and dissipation rate from running averages.
// Theoretical model only: no measured data. Pure functions and typed arrays, no DOM. Lattice units inside.
'use strict';

const LB_C = [[0, 0], [1, 0], [0, 1], [-1, 0], [0, -1], [1, 1], [-1, 1], [-1, -1], [1, -1]];
const LB_W = [4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36];
const LB_OPP = [0, 3, 4, 1, 2, 7, 8, 5, 6];

// Domain nx (along the flow, x) by ny (across, y). o.nu: kinematic viscosity; o.uIn(y): inflow velocity;
// o.walls: solid rows at y = 0 and y = ny - 1 (otherwise periodic in y).
function lbCreate(nx, ny, o) {
  const N = nx * ny;
  const s = {
    nx, ny, N, nu: o.nu, walls: !!o.walls, periodicX: !!o.periodicX, ramp: o.ramp || 0, cs: o.smagorinsky || 0, sponge: o.sponge || 0, uIn: o.uIn || (() => 0),
    f: new Float64Array(9 * N), g: new Float64Array(9 * N), solid: new Uint8Array(N),
    rho: new Float64Array(N), ux: new Float64Array(N), uy: new Float64Array(N), uxOld: new Float64Array(N), uyOld: new Float64Array(N),
    t: 0, stats: null,
  };
  const tau = 3 * o.nu + 0.5;
  s.tau0 = tau;
  s.wp = 1 / tau;
  s.wm = 1 / (0.5 + 0.25 / (tau - 0.5));
  if (s.walls) for (let x = 0; x < nx; x++) { s.solid[x] = 1; s.solid[(ny - 1) * nx + x] = 1; }
  for (let y = 0; y < ny; y++) for (let x = 0; x < nx; x++) {
    const k = y * nx + x;
    if (s.solid[k]) continue;
    lbSetEq(s, k, 1, o.u0 ? o.u0(x, y)[0] : s.ramp ? 0 : s.uIn(y), o.u0 ? o.u0(x, y)[1] : 0);
  }
  return s;
}
function lbSetEq(s, k, rho, ux, uy) {
  const uu = ux * ux + uy * uy;
  for (let i = 0; i < 9; i++) {
    const cu = LB_C[i][0] * ux + LB_C[i][1] * uy;
    s.f[9 * k + i] = LB_W[i] * rho * (1 + 3 * cu + 4.5 * cu * cu - 1.5 * uu);
  }
  s.rho[k] = rho; s.ux[k] = ux; s.uy[k] = uy;
}
// Mark a solid ellipse: centre (cx, cy), semi-axes ax (along x, the flow) and ay (across).
function lbAddEllipse(s, cx, cy, ax, ay) {
  for (let y = 0; y < s.ny; y++) for (let x = 0; x < s.nx; x++) {
    const dx = (x - cx) / ax, dy = (y - cy) / ay;
    if (dx * dx + dy * dy <= 1) { const k = y * s.nx + x; s.solid[k] = 1; s.ux[k] = 0; s.uy[k] = 0; for (let i = 0; i < 9; i++) s.f[9 * k + i] = 0; }
  }
  s.src = null;
}

// Neighbour table for streaming: src[9k+i] = index of the cell population i comes from, or -1 for a
// bounce-back from a solid, or -2 at the inlet and outlet ends (filled in afterwards). Rebuilt when solids change.
function lbLinks(s) {
  const { nx, ny, N, solid } = s, src = new Int32Array(9 * N);
  for (let y = 0; y < ny; y++) for (let x = 0; x < nx; x++) {
    const k = y * nx + x;
    for (let i = 0; i < 9; i++) {
      let sx = x - LB_C[i][0], sy = y - LB_C[i][1];
      if (s.periodicX) { if (sx < 0) sx += nx; else if (sx >= nx) sx -= nx; }
      if (sy < 0) sy += ny; else if (sy >= ny) sy -= ny;
      if (sx < 0 || sx >= nx) { src[9 * k + i] = -2; continue; }
      const ks = sy * nx + sx;
      src[9 * k + i] = solid[ks] ? -1 : 9 * ks + i;
    }
  }
  s.src = src; s.linkedSolids = solid.reduce((a, b) => a + b, 0);
}

// One time step: macroscopic values, collision, streaming with bounce-back, inlet and outlet.
function lbStep(s) {
  if (!s.src) lbLinks(s);
  const { nx, ny, N, solid, rho, ux, uy, src } = s;
  let wp = s.wp, wm = s.wm;
  const les = s.cs > 0 || s.sponge > 0, c36 = 36 * s.cs * s.cs, xs0 = nx - s.sponge;
  let tau0 = s.tau0;
  const f = s.f, g = s.g;
  s.uxOld.set(ux); s.uyOld.set(uy);
  const w0 = 4 / 9, w1 = 1 / 9, w2 = 1 / 36;
  for (let k = 0, x = 0; k < N; k++, x = x === nx - 1 ? 0 : x + 1) {
    if (solid[k]) continue;
    const b = 9 * k;
    const f0 = f[b], f1 = f[b + 1], f2 = f[b + 2], f3 = f[b + 3], f4 = f[b + 4], f5 = f[b + 5], f6 = f[b + 6], f7 = f[b + 7], f8 = f[b + 8];
    const r = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
    const u = (f1 - f3 + f5 - f6 - f7 + f8) / r, v = (f2 - f4 + f5 + f6 - f7 - f8) / r;
    rho[k] = r; ux[k] = u; uy[k] = v;
    const uu = 1.5 * (u * u + v * v);
    if (les) {
      // sponge layer before the outlet: the viscosity rises towards tau = 1 so that eddies leave without reflection
      tau0 = x > xs0 ? s.tau0 + (1 - s.tau0) * ((x - xs0) / s.sponge) ** 2 : s.tau0;
      // Smagorinsky eddy viscosity from the non-equilibrium momentum flux Q (Yu, Girimaji and Luo 2005):
      // tau = (tau0 + sqrt(tau0^2 + 4 C^2 cs^-4 |Q|)) / 2, |Q| = (Q_ab Q_ab)^(1/2), cs^-4 = 9
      const qxx = f1 + f3 + f5 + f6 + f7 + f8 - r * (1 / 3 + u * u), qyy = f2 + f4 + f5 + f6 + f7 + f8 - r * (1 / 3 + v * v), qxy = f5 - f6 + f7 - f8 - r * u * v;
      const tau = 0.5 * (tau0 + Math.sqrt(tau0 * tau0 + c36 * Math.sqrt(qxx * qxx + qyy * qyy + 2 * qxy * qxy)));
      wp = 1 / tau; wm = 1 / (0.5 + 0.25 / (tau - 0.5));
    }
    f[b] = f0 - wp * (f0 - w0 * r * (1 - uu));
    let cu, ep, em, dp, dm;
    cu = u; ep = w1 * r * (1 + 4.5 * cu * cu - uu); em = w1 * r * 3 * cu; dp = wp * (0.5 * (f1 + f3) - ep); dm = wm * (0.5 * (f1 - f3) - em); f[b + 1] = f1 - dp - dm; f[b + 3] = f3 - dp + dm;
    cu = v; ep = w1 * r * (1 + 4.5 * cu * cu - uu); em = w1 * r * 3 * cu; dp = wp * (0.5 * (f2 + f4) - ep); dm = wm * (0.5 * (f2 - f4) - em); f[b + 2] = f2 - dp - dm; f[b + 4] = f4 - dp + dm;
    cu = u + v; ep = w2 * r * (1 + 4.5 * cu * cu - uu); em = w2 * r * 3 * cu; dp = wp * (0.5 * (f5 + f7) - ep); dm = wm * (0.5 * (f5 - f7) - em); f[b + 5] = f5 - dp - dm; f[b + 7] = f7 - dp + dm;
    cu = v - u; ep = w2 * r * (1 + 4.5 * cu * cu - uu); em = w2 * r * 3 * cu; dp = wp * (0.5 * (f6 + f8) - ep); dm = wm * (0.5 * (f6 - f8) - em); f[b + 6] = f6 - dp - dm; f[b + 8] = f8 - dp + dm;
  }
  // streaming (pull) with halfway bounce-back
  for (let k = 0; k < N; k++) {
    if (solid[k]) continue;
    const b = 9 * k;
    for (let i = 0; i < 9; i++) {
      const q = src[b + i];
      g[b + i] = q >= 0 ? f[q] : q === -1 ? f[b + LB_OPP[i]] : f[b + i];
    }
  }
  // inlet
  // the inflow is raised smoothly over s.ramp steps, so the start sends no strong pressure waves
  const rampF = s.ramp ? 0.5 * (1 - Math.cos(Math.PI * Math.min(1, s.t / s.ramp))) : 1;
  for (let y = 0; y < ny && !s.periodicX; y++) {
    // inlet: prescribed velocity
    const k = y * nx;
    if (!solid[k]) {
      // Zou and He (1997): density from the known populations, then the three unknown ones (moving in +x)
      const u = rampF * s.uIn(y), q = 9 * k;
      const r = (g[q] + g[q + 2] + g[q + 4] + 2 * (g[q + 3] + g[q + 6] + g[q + 7])) / (1 - u);
      g[q + 1] = g[q + 3] + (2 / 3) * r * u;
      g[q + 5] = g[q + 7] - 0.5 * (g[q + 2] - g[q + 4]) + (1 / 6) * r * u;
      g[q + 8] = g[q + 6] + 0.5 * (g[q + 2] - g[q + 4]) + (1 / 6) * r * u;
    }
    // outlet: equilibrium at density 1 and the velocity of the cell upstream
    const ko = y * nx + nx - 1, kq = ko - 1;
    if (!solid[ko]) {
      const u = ux[kq], v = uy[kq], uu = 1.5 * (u * u + v * v);
      for (let i = 0; i < 9; i++) { const cu = LB_C[i][0] * u + LB_C[i][1] * v; g[9 * ko + i] = LB_W[i] * (1 + 3 * cu + 4.5 * cu * cu - uu); }
    }
  }
  s.f = g; s.g = f;
  s.t++;
}

// Velocity at a point by bilinear interpolation of a field (0 inside solids).
function lbInterp(s, F, x, y) {
  const { nx, ny } = s;
  if (x < 0) x = 0; if (x > nx - 1.001) x = nx - 1.001;
  const x0 = Math.floor(x), fx = x - x0;
  let y0 = Math.floor(y); const fy = y - y0;
  y0 = ((y0 % ny) + ny) % ny; const y1 = (y0 + 1) % ny;
  const a = F[y0 * nx + x0], b = F[y0 * nx + x0 + 1], c = F[y1 * nx + x0], d = F[y1 * nx + x0 + 1];
  return (1 - fy) * ((1 - fx) * a + fx * b) + fy * ((1 - fx) * c + fx * d);
}

// ---------- turbulence statistics ----------
// Running sums every call: u, v, u^2, v^2, uv, and the squared strain rate s_ij s_ij; the mean-flow strain
// is formed from the mean velocity at the end, so the dissipation rate of the fluctuations is
//   eps = 2 nu ( <s_ij s_ij> - S_ij S_ij ),   k = (<u'^2> + <v'^2>) / 2   (two components only).
function lbStatsReset(s) {
  const N = s.N;
  s.stats = { n: 0, u: new Float64Array(N), v: new Float64Array(N), uu: new Float64Array(N), vv: new Float64Array(N), uv: new Float64Array(N), ss: new Float64Array(N) };
}
function lbStrain(s, U, V, k) {
  // central differences, one-sided next to solids and at the ends
  const { nx, ny, solid } = s, x = k % nx, y = (k - x) / nx;
  const val = (F, xx, yy) => { yy = ((yy % ny) + ny) % ny; const q = yy * nx + xx; return solid[q] ? 0 : F[q]; };
  const dx = (F) => (x > 0 && x < nx - 1 ? (val(F, x + 1, y) - val(F, x - 1, y)) / 2 : 0);
  const dy = (F) => (val(F, x, y + 1) - val(F, x, y - 1)) / 2;
  const ux = dx(U), uy = dy(U), vx = dx(V), vy = dy(V);
  return ux * ux + vy * vy + 0.5 * (uy + vx) * (uy + vx);
}
function lbStatsAdd(s) {
  if (!s.stats) lbStatsReset(s);
  const st = s.stats, { N, ux, uy, solid } = s;
  for (let k = 0; k < N; k++) {
    if (solid[k]) continue;
    const u = ux[k], v = uy[k];
    st.u[k] += u; st.v[k] += v; st.uu[k] += u * u; st.vv[k] += v * v; st.uv[k] += u * v;
    st.ss[k] += lbStrain(s, ux, uy, k);
  }
  st.n++;
}
// Mean speed, turbulent kinetic energy and dissipation rate per cell.
function lbStatsFields(s) {
  const st = s.stats, N = s.N, n = Math.max(1, st.n);
  const U = new Float64Array(N), V = new Float64Array(N), K = new Float64Array(N), E = new Float64Array(N);
  for (let k = 0; k < N; k++) { U[k] = st.u[k] / n; V[k] = st.v[k] / n; K[k] = Math.max(0, 0.5 * (st.uu[k] / n - U[k] * U[k] + st.vv[k] / n - V[k] * V[k])); }
  for (let k = 0; k < N; k++) if (!s.solid[k]) E[k] = Math.max(0, 2 * s.nu * (st.ss[k] / n - lbStrain(s, U, V, k)));
  return { U, V, K, E, n: st.n };
}

// ---------- inertial particles ----------
// dv/dt = beta Du/Dt + (u - v) / tp,  beta = 3 / (2 rho_p/rho + 1), dx/dt = v   (lattice units, time step 1)
// Du/Dt from the change of the field over the last step and the convective term (u . grad) u.
function lbMoveParticle(s, p, beta, tp) {
  const u = lbInterp(s, s.ux, p.x, p.y), v = lbInterp(s, s.uy, p.x, p.y);
  const uo = lbInterp(s, s.uxOld, p.x, p.y), vo = lbInterp(s, s.uyOld, p.x, p.y);
  const h = 0.5;
  const dudx = (lbInterp(s, s.ux, p.x + h, p.y) - lbInterp(s, s.ux, p.x - h, p.y)) / (2 * h), dudy = (lbInterp(s, s.ux, p.x, p.y + h) - lbInterp(s, s.ux, p.x, p.y - h)) / (2 * h);
  const dvdx = (lbInterp(s, s.uy, p.x + h, p.y) - lbInterp(s, s.uy, p.x - h, p.y)) / (2 * h), dvdy = (lbInterp(s, s.uy, p.x, p.y + h) - lbInterp(s, s.uy, p.x, p.y - h)) / (2 * h);
  const ax = (u - uo) + u * dudx + v * dudy, ay = (v - vo) + u * dvdx + v * dvdy;
  // drag treated exactly over the step for small tp (exponential integrator)
  const e = Math.exp(-1 / tp);
  const vx = u + (p.vx - u) * e + beta * ax * tp * (1 - e), vy = v + (p.vy - v) * e + beta * ay * tp * (1 - e);
  p.x += 0.5 * (p.vx + vx); p.y += 0.5 * (p.vy + vy); p.vx = vx; p.vy = vy;
  if (p.y < 0) p.y += s.ny; else if (p.y >= s.ny) p.y -= s.ny;
}

if (typeof module !== 'undefined') {
  module.exports = { LB_C, LB_W, LB_OPP, lbLinks, lbCreate, lbSetEq, lbAddEllipse, lbStep, lbInterp, lbStatsReset, lbStatsAdd, lbStatsFields, lbStrain, lbMoveParticle };
}
