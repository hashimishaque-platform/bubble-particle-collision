# Tutorial 1: collision efficiency of quartz particles

**Question.** How likely is a quartz particle of 10, 30 or 60 µm to hit a bubble of 0.6 or
1.0 mm in water with frother, and how much does particle inertia matter?

## In the page

1. Set the bubble diameter to 1.0 mm, the surface to "with frother (immobile)", the material
   to quartz and the diameter to 30 µm. Press Compute.
2. Figure 1 shows the streamlines and the paths of particles; the thick red path is the
   grazing one. Figure 2 shows the collision efficiency E across particle sizes.
3. Press "Results (JSON)" and "Figure 2 curves (CSV)" to keep the numbers, and "Link to these
   settings" to keep the case.

## From the command line

```
node cli.js examples/collision-quartz.json --csv quartz.csv
```

This takes about five seconds. Part of the result:

| Bubble (mm) | Re | Particle (µm) | Stokes number | E, interception | E, with settling | E, with inertia | Weber-Paddock |
|---|---|---|---|---|---|---|---|
| 0.6 | 40.2 | 30 | 0.035 | 0.0125 | 0.0241 | 0.0242 | 0.0133 |
| 1.0 | 111.8 | 10 | 0.0039 | 0.00083 | 0.00153 | 0.00154 | 0.00085 |
| 1.0 | 111.8 | 30 | 0.035 | 0.00728 | 0.0138 | 0.0135 | 0.00765 |
| 1.0 | 111.8 | 60 | 0.135 | 0.0281 | 0.0540 | 0.0570 | 0.0306 |

## How to read it

- Interception alone agrees with the Weber and Paddock correlation within about 8 %, as the
  checks show for Re = 1 to 100.
- Settling nearly doubles E for quartz, because the particles sink across the streamlines
  while the bubble rises through them.
- Inertia changes little at a Stokes number of 0.035 but raises E by about 6 % at 0.135.
- E falls with bubble size for fine particles: the smaller bubble collects better.

The value 0.01346 for 30 µm quartz and a 1.0 mm bubble is one of the program's checks: the
Runge-Kutta and exponential integrators give it to 0.02 %.
