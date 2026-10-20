#!/usr/bin/env node
// Batch runs without the browser: node cli.js case.json [--out result.json] [--csv result.csv]
// The same code as the page (flow.js, swarm.js, lbm.js). See examples/*.json and docs/cli.md.
'use strict';
const fs = require('fs');
const F = require('./flow.js'), S = require('./swarm.js'), LB = require('./lbm.js');
const VERSION = (fs.readFileSync(__dirname + '/build.py', 'utf8').match(/VERSION = '([^']+)'/) || [])[1] || 'unknown';
const RHO_L = 998.2, MU = 1.002e-3, RHO_G = 1.2;
const list = (v) => (Array.isArray(v) ? v : [v]);

// collision efficiency for each bubble size and particle size: the computation behind Figures 1 and 2
function runCollision(c) {
  const surf = c.surface === 'mobile' ? 'free' : 'rigid', rows = [];
  for (const dbMm of list(c.bubbleDiameter_mm)) {
    const db = dbMm * 1e-3;
    if (db > (surf === 'free' ? 0.7e-3 : 1.2e-3)) throw new Error(`bubble diameter ${dbMm} mm is outside the range of the flow solver`);
    const U = F.terminalVelocity(db, RHO_L - RHO_G, RHO_L, MU, surf === 'rigid' ? F.cdRigid : F.cdClean), Re = RHO_L * U * db / MU;
    const st = F.solveFlow(Re, surf, Object.assign({ tol: c.tolerance || 1e-8, maxIt: 100000 }, F.gridFor(Re, surf, c.grid || 'normal')));
    const ev = F.psiEvaluator(st), sample = F.flowSampler(st);
    for (const dpUm of list(c.particleDiameter_um)) for (const rhop of list(c.particleDensity_kg_m3)) {
      const pp = F.particleParams({ db, U, rhoL: RHO_L, mu: MU }, { dp: dpUm * 1e-6, rhop });
      const Ei = F.collisionStream(ev, pp.rp, 0).E, Eg = F.collisionStream(ev, pp.rp, pp.ws).E;
      const Einr = c.inertia !== false && pp.St >= 1e-3 ? F.collisionInertial(sample, ev, pp).E : null;
      rows.push({ bubbleDiameter_mm: dbMm, surface: c.surface || 'immobile', riseSpeed_m_s: U, Re, dragCoefficient: F.dragCoefficient(st).cd,
        particleDiameter_um: dpUm, particleDensity_kg_m3: rhop, StokesNumber: pp.St, E_interception: Ei, E_interception_settling: Eg, E_withInertia: Einr,
        E_correlation: surf === 'rigid' ? F.ecWeberPaddock(pp.rp, Re) : F.ecSutherland(pp.rp), correlation: surf === 'rigid' ? 'Weber and Paddock (1983)' : 'Sutherland (1948)' });
    }
  }
  return rows;
}
// rise speed and shape from the correlations: Figure 6
function runRise(c) {
  const rows = [];
  for (const dMm of list(c.bubbleDiameter_mm)) for (const sig of list(c.surfaceTension_mN_m || 72.75)) {
    const d = dMm * 1e-3, s = sig * 1e-3, Eo = S.eotvos(d, s);
    rows.push({ bubbleDiameter_mm: dMm, surfaceTension_mN_m: sig, Eo, U_pure_m_s: S.riseVelocity(d, 'pure', S.SW.mu, s), U_slight_m_s: S.riseVelocity(d, 'slight', S.SW.mu, s), U_full_m_s: S.riseVelocity(d, 'full', S.SW.mu, s), aspectRatio_Wellek: S.aspectWellek(Eo) });
  }
  return rows;
}
// the live two-dimensional simulation of Figures 7 and 8, run for a number of steps
function runLive(c) {
  const D = 20, nx = 300, ny = 120, Ul = 0.04, xb = 180, sponge = 40;
  const d = c.bubbleDiameter_mm * 1e-3, U = c.speed_m_s, Re = Math.min(U * d / (MU / RHO_L), 1500), nu = Ul * D / Re, E = c.aspectRatio || 0.5;
  const s = LB.lbCreate(nx, ny, { nu, ramp: 1500, smagorinsky: 0.15, sponge, uIn: () => Ul });
  const R = D / 2, cx = R * Math.sqrt(E), ay = R / Math.sqrt(E), yb = ny / 2 + 0.3;
  LB.lbAddEllipse(s, xb, yb, cx, ay);
  const rods = [];
  if (c.gridDiametersUpstream) { const M = 1.25 * D, dr = D / 3, xg = xb - c.gridDiametersUpstream * D; for (let y = M / 2; y < ny; y += M) { LB.lbAddEllipse(s, xg, y, dr / 2, dr / 2); rods.push([xg, y, dr / 2]); } }
  const dpl = (c.particleDiameter_um || 33) * 1e-6 / d * D, ratio = (c.particleDensity_kg_m3 || 1050) / RHO_L;
  const tp = Math.max(1e-3, (ratio + 0.5) * dpl * dpl / (18 * nu)), beta = 3 / (2 * ratio + 1);
  const rand = S.swRandom(c.seed || 1), steps = c.steps || 20000, warm = c.warmup || 5000, hits = [];
  let parts = [];
  const x0 = rods.length ? Math.min(rods[0][0] + D, xb - 3 * D) : xb - 3 * D;
  for (let n = 0; n < steps; n++) {
    LB.lbStep(s);
    if (s.t > warm && s.t % 4 === 0) LB.lbStatsAdd(s);
    if (s.t < 1500) continue;
    const y = yb + (2 * rand() - 1) * 1.3 * ay, x = x0 + rand();
    parts.push({ x, y, vx: LB.lbInterp(s, s.ux, x, y), vy: LB.lbInterp(s, s.uy, x, y) });
    const ex = cx + dpl / 2, ey = ay + dpl / 2, keep = [];
    for (const p of parts) {
      LB.lbMoveParticle(s, p, beta, tp);
      const dx = (p.x - xb) / ex, dy = (p.y - yb) / ey;
      if (dx * dx + dy * dy <= 1) { hits.push(Math.atan2(Math.abs(dy), -dx) * 180 / Math.PI); continue; }
      if (rods.some((r) => (p.x - r[0]) ** 2 + (p.y - r[1]) ** 2 <= (r[2] + dpl / 2) ** 2)) continue;
      if (p.x < nx - sponge) keep.push(p);
    }
    parts = keep;
  }
  const st = s.stats && s.stats.n > 2 ? LB.lbStatsFields(s) : null, kS = U * U / (Ul * Ul), eS = (U ** 3 / d) / (Ul ** 3 / D), ring = [];
  if (st) for (let a = 0; a <= 180; a += 5) {
    const t = a * Math.PI / 180, x = xb - 1.5 * cx * Math.cos(t), dy = 1.5 * ay * Math.sin(t);
    const v = (Fd) => 0.5 * (LB.lbInterp(s, Fd, x, yb + dy) + LB.lbInterp(s, Fd, x, yb - dy));
    ring.push({ angle_deg: a, k_m2_s2: v(st.K) * kS, eps_W_kg: v(st.E) * eS });
  }
  return [{ Re, StokesNumber: tp * Ul / D, steps, averages: st ? st.n : 0, collisionsFront: hits.filter((q) => q < 90).length, collisionsRear: hits.filter((q) => q >= 90).length, collisionAngles_deg: hits, ring }];
}

function toCSV(rows) {
  const flat = rows.map((r) => Object.fromEntries(Object.entries(r).filter(([, v]) => !Array.isArray(v) && typeof v !== 'object' || v === null)));
  const head = [...new Set(flat.flatMap((r) => Object.keys(r)))];
  return head.join(',') + '\n' + flat.map((r) => head.map((h) => (r[h] === null || r[h] === undefined ? '' : typeof r[h] === 'number' ? String(+r[h].toPrecision(10)) : `"${String(r[h]).replace(/"/g, '""')}"`)).join(',')).join('\n') + '\n';
}

if (require.main === module) {
  const args = process.argv.slice(2);
  if (!args.length || args[0] === '--help') {
    console.log('Usage: node cli.js case.json [--out result.json] [--csv result.csv]\nTools: "collision", "rise", "live". Examples in examples/.');
    process.exit(args.length ? 0 : 1);
  }
  const c = JSON.parse(fs.readFileSync(args[0], 'utf8'));
  const run = { collision: runCollision, rise: runRise, live: runLive }[c.tool];
  if (!run) { console.error(`unknown tool "${c.tool}"; use collision, rise or live`); process.exit(1); }
  const t0 = Date.now(), rows = run(c);
  const out = { program: 'Bubble and particle collision', version: VERSION, case: c, seconds: (Date.now() - t0) / 1000, results: rows };
  const o = args.indexOf('--out'), v = args.indexOf('--csv');
  if (o >= 0) fs.writeFileSync(args[o + 1], JSON.stringify(out, null, 1));
  if (v >= 0) fs.writeFileSync(args[v + 1], toCSV(rows));
  if (o < 0 && v < 0) console.log(JSON.stringify(out, null, 1));
}
module.exports = { runCollision, runRise, runLive, toCSV };
