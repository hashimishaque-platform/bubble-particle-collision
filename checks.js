// Runs every check shown on the page and writes checks.json. Usage: node checks.js
// Theoretical checks only: the expected values are exact solutions or published equations and correlations.
'use strict';
const F = require('./flow.js'), P = require('./piv.js');
const fs = require('fs');
const rows = [];
const add = (group, label, expected, normal, fine, source) => {
  rows.push({ group, label, expected, normal, fine, source });
  console.log(group, '|', label, '|', expected, '|', normal, '|', fine);
};
const fx = (x, d) => Number(x).toFixed(d);
const sup = (s) => `<sup>${s}</sup>`;
const sci = (x, d) => { const e = Math.floor(Math.log10(Math.abs(x))); return `${(x / 10 ** e).toFixed(d)}\u00a0×\u00a010${sup(String(e).replace('-', '−'))}`; };
const pct = (a, b) => { const v = 100 * (a / b - 1); return `${v < 0 ? '−' : '+'}${Math.abs(v).toFixed(2)}\u00a0%`; };
const RATIO = 'd<sub>p</sub>/d<sub>b</sub>';
const cache = new Map();
const solve = (Re, surf, q) => {
  const key = `${Re}|${surf}|${q}`;
  if (!cache.has(key)) cache.set(key, F.solveFlow(Re, surf, Object.assign({ tol: 1e-9, maxIt: 100000 }, F.gridFor(Re, surf, q))));
  return cache.get(key);
};

// ---- flow, immobile surface ----
const G1 = 'Flow past a sphere with an immobile surface';
for (const Re of [0.1, 50, 100, 150]) {
  const a = solve(Re, 'rigid', 'normal'), b = solve(Re, 'rigid', 'fine');
  const ca = F.dragCoefficient(a).cd, cb = F.dragCoefficient(b).cd;
  if (Re === 0.1) {
    const pp = 24 / Re * (1 + 3 * Re / 16 + 9 / 160 * Re * Re * Math.log(Re / 2));
    add(G1, `Drag coefficient, Re = ${Re}`, fx(pp, 1), fx(ca, 1), fx(cb, 1), 'Proudman and Pearson (1957) expansion');
    continue;
  }
  const jp = { 50: [1.57, 0.40], 100: [1.08, 0.88], 150: [0.90, 1.20] }[Re];
  add(G1, `Drag coefficient, Re = ${Re}`, fx(jp[0], 2), fx(ca, 3), fx(cb, 3), 'Johnson and Patel (1999)');
  add(G1, `Wake length / diameter, Re = ${Re}`, fx(jp[1], 2), fx(F.wakeLength(a), 2), fx(F.wakeLength(b), 2), 'Johnson and Patel (1999)');
  add(G1, `Separation angle from the rear, Re = ${Re}`, fx(180 - F.sepAngleClift(Re), 1) + '°', fx(180 - F.separationAngle(a), 1) + '°', fx(180 - F.separationAngle(b), 1) + '°',
    Re === 100 ? 'Clift, Grace and Weber (1978); Johnson and Patel (1999): 53°' : 'Clift, Grace and Weber (1978)');
}

// ---- flow, mobile surface ----
const G2 = 'Flow past a spherical bubble with a clean (mobile) surface';
for (const Re of [0.1, 10, 50, 100]) {
  const a = solve(Re, 'free', 'normal'), b = solve(Re, 'free', 'fine');
  const d = Re < 1 ? 1 : 3;
  add(G2, `Drag coefficient, Re = ${Re}`, fx(F.cdClean(Re), d), fx(F.dragCoefficient(a).cd, d), fx(F.dragCoefficient(b).cd, d),
    Re < 1 ? 'Mei, Klausner and Lawrence (1994); Hadamard-Rybczynski limit 16/Re = 160' : 'Mei, Klausner and Lawrence (1994)');
}

// ---- collision, exact limits (analytic stream functions) ----
const G3 = 'Collision efficiency, exact limits';
{
  const rp = 0.05, r = 1 + rp;
  add(G3, `Stokes flow, interception, ${RATIO} = 0.05`, sci(r * r - 1.5 * r + 0.5 / r, 4), sci(F.collisionStream(F.ANALYTIC.stokesRigid, rp, 0).E, 4), '', 'Exact for Stokes flow; Gaudin (1957) for small particles');
  add(G3, `Potential flow, interception, ${RATIO} = 0.05`, fx(r * r - 1 / r, 5), fx(F.collisionStream(F.ANALYTIC.potential, rp, 0).E, 5), '', 'Exact for potential flow; Sutherland (1948) for small particles');
  const ws = 0.1;
  add(G3, 'Stokes flow, settling only (point particle, settling speed 0.1 U)', fx(ws / (1 + ws), 5), fx(F.collisionStream(F.ANALYTIC.stokesRigid, 1e-6, ws).E, 5), '', 'Exact: w<sub>s</sub> / (U + w<sub>s</sub>)');
}

// ---- collision, computed flows ----
const G4 = 'Collision efficiency, immobile surface, interception only';
for (const Re of [1, 10, 50, 100]) {
  const a = solve(Re, 'rigid', 'normal'), b = solve(Re, 'rigid', 'fine');
  const rp = 0.005;
  const ea = F.collisionStream(F.psiEvaluator(a), rp, 0).E / (rp * rp), eb = F.collisionStream(F.psiEvaluator(b), rp, 0).E / (rp * rp);
  const wp = 1.5 * (1 + (3 / 16) * Re / (1 + 0.249 * Math.pow(Re, 0.56)));
  const yl = 1.5 + 4 * Math.pow(Re, 0.72) / 15;
  add(G4, `E / (${RATIO})² for small particles, Re = ${Re}`, `${fx(wp, 2)} (Yoon and Luttrell: ${fx(yl, 2)})`, fx(ea, 2), fx(eb, 2), 'Weber and Paddock (1983)');
}

// ---- particle trajectories against the stream function ----
const G5 = 'Particle trajectories';
for (const surf of ['rigid', 'free']) {
  const out = [];
  for (const q of ['normal', 'fine']) {
    const st = solve(100, surf, q), ev = F.psiEvaluator(st), sample = F.flowSampler(st);
    const pp = { rp: 0.03, ws: 0.02, St: 1e-7, beta: 0 };
    out.push([F.collisionStream(ev, pp.rp, pp.ws).E, F.collisionInertial(sample, ev, pp).E]);
  }
  add(G5, `Without inertia (Stokes number 10${sup('−7')}), ${surf === 'rigid' ? 'immobile' : 'mobile'} surface, Re = 100, ${RATIO} = 0.03, w<sub>s</sub>/U = 0.02: E from integrated paths`,
    `${fx(out[0][0], 5)} (fine grid: ${fx(out[1][0], 5)})`, fx(out[0][1], 5), fx(out[1][1], 5),
    `Exact for the same flow, from the particle stream function; difference ${pct(out[0][1], out[0][0])}`);
}
{
  // 30 um quartz particle and a 1.0 mm bubble with frother, as on the page: both integrators
  const RHO_L = 998.2, MU = 1.002e-3, db = 1e-3;
  const U = F.terminalVelocity(db, RHO_L - 1.2, RHO_L, MU, F.cdRigid), Re = RHO_L * U * db / MU;
  const st = F.solveFlow(Re, 'rigid', Object.assign({ tol: 1e-9, maxIt: 100000 }, F.gridFor(Re, 'rigid', 'normal')));
  const ev = F.psiEvaluator(st), sample = F.flowSampler(st);
  const pp = F.particleParams({ db, U, rhoL: RHO_L, mu: MU }, { dp: 30e-6, rhop: 2650 });
  const a = F.collisionInertial(sample, ev, pp, { method: 'rk4' }), b = F.collisionInertial(sample, ev, pp, { method: 'exp' });
  add(G5, `With inertia: 30 µm quartz, 1.0 mm bubble with frother (Re = ${fx(Re, 1)}, Stokes number ${fx(pp.St, 3)}), exponential integrator against Runge-Kutta`,
    fx(a.E, 5), fx(b.E, 5), '', `Same flow; difference ${pct(b.E, a.E)}`);
}

// ---- PIV ----
{
  const w = 256, h = 256, rand = P.mulberry32(7), dx = 2.3, dy = -1.7;
  const nP = Math.round(0.03 * (w + 40) * (h + 40)), A = [], B = [];
  for (let k = 0; k < nP; k++) {
    const x = -20 + rand() * (w + 40), y = -20 + rand() * (h + 40), I0 = 60 + 140 * rand();
    A.push([x - dx / 2, y - dy / 2, I0]); B.push([x + dx / 2, y + dy / 2, I0]);
  }
  const ia = P.finishImage(P.renderImage(w, h, A, 3), w, h, rand), ib = P.finishImage(P.renderImage(w, h, B, 3), w, h, rand);
  const r = P.evaluate(ia, ib, w, h, {});
  let s = 0, n = 0;
  for (let k = 0; k < r.nx * r.ny; k++) { if (!r.valid[k]) continue; s += (r.U[k] - dx) ** 2 + (r.V[k] - dy) ** 2; n++; }
  add('Synthetic PIV', 'Uniform shift of 2.3 and −1.7 pixels, 32 px windows: random error of one displacement component', 'about 0.05 px or less', fx(Math.sqrt(s / n / 2), 3) + ' px', '', 'Typical for 32 px windows: Raffel et al. (2018)');
}

fs.writeFileSync('checks.json', JSON.stringify(rows, null, 1));
console.log('written', rows.length, 'rows');
