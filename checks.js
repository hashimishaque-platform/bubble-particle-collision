// Runs every check shown on the page and writes checks.json. Usage: node checks.js
// Theoretical checks only: the expected values are exact solutions or published equations and correlations.
'use strict';
const F = require('./flow.js'), P = require('./piv.js'), S = require('./swarm.js');
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

// ---- bubble swarm and turbulence ----
{
  const G = 'Bubble swarm and synthetic turbulence';
  const rand = S.swRandom(11);
  // a spherical bubble of radius 1 mm moving with (0.05, 0.3, -0.02) m/s
  const b = { x: 0, y: 0, z: 0, vx: 0.05, vy: 0.3, vz: -0.02, a: 1, c: 1, alpha: 0 };
  let eBC = 0, eDiv = 0, gMax = 0;
  const h = 1e-5;
  for (let i = 0; i < 400; i++) {
    const zz = 2 * rand() - 1, ph = 2 * Math.PI * rand(), sq = Math.sqrt(1 - zz * zz), n = [sq * Math.cos(ph), sq * Math.sin(ph), zz];
    const r = 1 + 1e-9, u = S.bubbleFlow([b], r * n[0], r * n[1], r * n[2]);
    const un = u[0] * n[0] + u[1] * n[1] + u[2] * n[2], Un = b.vx * n[0] + b.vy * n[1] + b.vz * n[2];
    eBC = Math.max(eBC, Math.abs(un - Un) / Math.hypot(b.vx, b.vy, b.vz));
    const R = 1.2 + 3 * rand(), x = R * n[0], y = R * n[1], z = R * n[2], f = (p, q, w) => S.bubbleFlow([b], p, q, w);
    const dudx = (f(x + h, y, z)[0] - f(x - h, y, z)[0]) / (2 * h);
    const div = dudx + (f(x, y + h, z)[1] - f(x, y - h, z)[1]) / (2 * h) + (f(x, y, z + h)[2] - f(x, y, z - h)[2]) / (2 * h);
    eDiv = Math.max(eDiv, Math.abs(div)); gMax = Math.max(gMax, Math.abs(dudx));
  }
  add(G, 'Flow around one bubble: velocity normal to the surface minus the bubble velocity normal to it, largest of 400 points, relative to the bubble speed', '0', sci(eBC + 1e-16, 1), '', 'Kinematic condition at the surface, exact for the dipole field (Batchelor 1967)');
  add(G, 'Flow around one bubble: divergence, largest of 400 points, relative to the largest velocity gradient', '0 (round-off)', sci(eDiv / gMax, 1), '', 'Incompressible liquid');
  // turbulence at 1 W/kg, L = 20 mm, eddies down to 2 mm
  const T = S.makeTurbulence({ eps: 1, L: 0.02, kmax: 1 / 0.2e-3, modes: 200, seed: 5 });
  let e = 0, dMax = 0, gM = 0;
  const N = 40000, hh = 1e-7;
  for (let i = 0; i < N; i++) { const u = S.turbulentFlow(T, rand() * 0.4, rand() * 0.4, rand() * 0.4); e += 0.5 * (u[0] ** 2 + u[1] ** 2 + u[2] ** 2); }
  for (let i = 0; i < 100; i++) {
    const x = rand() * 0.02, y = rand() * 0.02, z = rand() * 0.02, f = (p, q, w) => S.turbulentFlow(T, p, q, w);
    const dudx = (f(x + hh, y, z)[0] - f(x - hh, y, z)[0]) / (2 * hh);
    const div = dudx + (f(x, y + hh, z)[1] - f(x, y - hh, z)[1]) / (2 * hh) + (f(x, y, z + hh)[2] - f(x, y, z - hh)[2]) / (2 * hh);
    dMax = Math.max(dMax, Math.abs(div)); gM = Math.max(gM, Math.abs(dudx));
  }
  add(G, 'Synthetic turbulence, 1 W/kg, L = 20 mm, 200 modes: mean kinetic energy over 40 000 points divided by the sum of E(k) dk over the modes', '1 (sampling error about 1 %)', fx(e / N / T.energy, 4), '', 'Kinematic simulation, Fung et al. (1992)');
  add(G, 'Same: sum of E(k) dk over the modes divided by the integral of the Kolmogorov spectrum over the same wavenumbers', '1', fx(T.energy / T.exact, 6), '', 'Band sum against the exact integral');
  add(G, 'Same: divergence, largest of 100 points, relative to the largest velocity gradient', '0 (round-off)', sci(dMax / gM, 1), '', 'Each mode is perpendicular to its wave vector');
  // rise speed: small bubble against Schiller-Naumann, large bubble against the shape branch in closed form
  const d1 = 0.1e-3, u1 = S.riseVelocity(d1, 'full'), u2 = F.terminalVelocity(d1, S.SW.rhoL - S.SW.rhoG, S.SW.rhoL, S.SW.mu, F.cdRigid);
  add(G, 'Rise speed of a 0.1 mm bubble, fully contaminated (Tomiyama et al. 1998) divided by Schiller and Naumann (1933)', '1', fx(u1 / u2, 6), '', 'Tomiyama et al. reduce to Schiller and Naumann at low Re and Eo');
  const d2 = 10e-3, Eo = S.eotvos(d2), cd = 8 / 3 * Eo / (Eo + 4), u3 = Math.sqrt(4 * (S.SW.rhoL - S.SW.rhoG) * S.SW.g * d2 / (3 * cd * S.SW.rhoL));
  add(G, 'Rise speed of a 10 mm bubble (m/s), slightly contaminated', `${fx(u3, 4)} (C<sub>D</sub> = 8 Eo / (3 (Eo + 4)), Eo = ${fx(Eo, 2)})`, fx(S.riseVelocity(d2, 'slight'), 4), '', 'Shape branch of Tomiyama et al. (1998), closed form');
  // lower surface tension: a 2.3 mm bubble at 29.7 mN/m, fully contaminated, lies on the shape branch
  const sg = 0.0297, d3 = 2.3e-3, Eo3 = S.eotvos(d3, sg), cd3 = 8 / 3 * Eo3 / (Eo3 + 4), u4 = Math.sqrt(4 * (S.SW.rhoL - S.SW.rhoG) * S.SW.g * d3 / (3 * cd3 * S.SW.rhoL));
  const u5 = S.riseVelocity(d3, 'full', S.SW.mu, sg);
  add(G, 'Rise speed of a 2.3 mm bubble at a surface tension of 29.7 mN/m, fully contaminated (m/s)', `${fx(u4, 4)} (shape branch, Eo = ${fx(Eo3, 2)})`, fx(u5, 4), '', 'Tomiyama et al. (1998), closed form');
  // tracks file reader and the time-step estimate, on a chain made and moved by the program
  const ch = S.makeSwarm({ n: 4, d: 2.5e-3, sys: 'full', chain: true, spacing: 10, box: [30, 40, 10] });
  let csv = 'LfdNr,TimeStep,ID,X,Y,Z,VX,VY,VZ,a,c,alpha\n', row = 0;
  for (let k = 1; k <= 60; k++) for (const q of S.advanceSwarm(ch.bubbles, (k - 1) * 0.5e-3, 1e6)) csv += [++row, k, q.id, q.x, q.y, q.z, q.vx, q.vy, q.vz, q.a, q.c, q.alpha].join(',') + '\n';
  const tr = S.tracksById(S.parseTracks(csv)), dtEst = S.trackTimeStep(tr), st = S.trackStats(tr.get(1));
  add(G, 'Tracks made by the program, 0.5 ms per time step: time step estimated from distance moved and velocity (ms)', '0.5', fx(dtEst * 1e3, 6), '', 'Exact for bubbles moving at constant speed');
  add(G, 'Same tracks: mean rise speed and diameter of one track', `${fx(ch.U, 4)} m/s, 2.500 mm`, `${fx(st.U, 4)} m/s, ${fx(st.d * 1e3, 3)} mm`, '', 'Values used to make the tracks');
  if (Math.abs(u5 / u4 - 1) > 1e-6 || Math.abs(dtEst / 0.5e-3 - 1) > 1e-7 || Math.abs(st.U / ch.U - 1) > 1e-9) { console.error('swarm track checks failed'); process.exitCode = 1; }
  if (eBC > 1e-6 || eDiv / gMax > 1e-4 || Math.abs(e / N / T.energy - 1) > 0.03 || dMax / gM > 1e-4 || Math.abs(u1 / u2 - 1) > 1e-6) { console.error('swarm checks failed'); process.exitCode = 1; }
}

fs.writeFileSync('checks.json', JSON.stringify(rows, null, 1));
console.log('written', rows.length, 'rows');
