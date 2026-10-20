# Tutorial 2: turbulence behind a bubble in a stream

**Question.** Where around a 2.3 mm bubble is the turbulence strongest, and where do fine
particles hit it, when the oncoming flow has passed a grid?

## In the page

1. In "Turbulence around a bubble: live flow simulation", keep 2.3 mm, 0.18 m/s (Re = 412),
   aspect ratio 0.5 and the grid 6 diameters upstream. Particles: 33 µm, 1050 kg/m³.
2. Press Play and set the speed to "fast". The rods shed rows of small eddies; behind the
   bubble larger eddies are shed in turn.
3. After about 5000 time steps the averages start. Switch the colour to "turbulent kinetic
   energy, averaged" and then to "dissipation rate, averaged".
4. Figure 8 shows k and ε around the bubble and where the particles hit.

## From the command line

```
node cli.js examples/live-grid.json --out live.json
```

20 000 steps took 84 s. The result had 3750 averages and 217 collisions:

| Angle from the front (°) | 0 | 30 | 60 | 90 | 120 | 150 | 180 |
|---|---|---|---|---|---|---|---|
| k (m²/s²) | 0.0015 | 0.0023 | 0.0050 | 0.0070 | 0.0135 | 0.0093 | 0.0205 |
| ε (W/kg) | 0.013 | 0.021 | 0.023 | 0.025 | 0.080 | 0.165 | 0.132 |

| Collision angle (°) | 0 to 20 | 20 to 40 | 40 to 60 | 60 to 80 | 80 to 100 | 100 to 120 | 120 to 140 | 140 to 160 | 160 to 180 |
|---|---|---|---|---|---|---|---|---|---|
| Collisions | 24 | 20 | 41 | 76 | 15 | 1 | 5 | 18 | 17 |

## How to read it

- The turbulent kinetic energy is more than ten times higher behind the bubble than in front of
  it, and the dissipation rate rises sharply past 120°. The wake, not the oncoming
  turbulence, dominates near the bubble.
- 176 of 217 collisions (81 %) are on the front half, most of them between 40° and 80°, where
  streamlines pass close to the surface. 41 (19 %) are on the rear half, between 140° and
  180°, where the shed eddies bring particles back to the bubble.
- Run it again with `"gridDiametersUpstream": 0` to see how much of the rear share comes from
  the oncoming turbulence.

The simulation is two-dimensional and the surface is immobile (see "Range of validity" on the
page), so read the numbers as an illustration of the mechanism.
