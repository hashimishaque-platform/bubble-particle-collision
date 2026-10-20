# Command-line use

`cli.js` runs the same code as the page from a case file, without a browser, so parameter
studies can be scripted. It needs Node.js 18 or newer and nothing else.

```
node cli.js case.json                         # results as JSON on the screen
node cli.js case.json --out result.json       # JSON file: version, case and results
node cli.js case.json --csv result.csv        # one row per result
```

Every output records the program version and the full case, so a result can be traced back.
Values given as lists are combined: every bubble size with every particle size and density.

## Tools

### `collision`: collision efficiency (Figures 1 and 2)

| Key | Meaning |
|---|---|
| `bubbleDiameter_mm` | number or list; up to 1.2 for an immobile surface, 0.7 for a mobile one |
| `surface` | `"immobile"` (with frother) or `"mobile"` (clean water) |
| `grid` | `"normal"` or `"fine"` |
| `particleDiameter_um`, `particleDensity_kg_m3` | numbers or lists |
| `inertia` | `true` (default) to trace paths with inertia when the Stokes number is at least 0.001 |

Output per case: rise speed, Re, drag coefficient, Stokes number, the collision efficiency by
interception, with settling, with inertia, and the correlation for comparison (Weber and
Paddock for an immobile surface, Sutherland for a mobile one).

### `rise`: rise speed and shape (Figure 6)

`bubbleDiameter_mm` and `surfaceTension_mN_m` (numbers or lists). Output: Eötvös number, rise
speed after Tomiyama et al. (1998) for pure, slightly and fully contaminated water, and the
aspect ratio after Wellek et al. (1966).

### `live`: two-dimensional flow past a bubble (Figures 7 and 8)

| Key | Meaning |
|---|---|
| `bubbleDiameter_mm`, `speed_m_s` | set the Reynolds number in water at 20 °C |
| `aspectRatio` | 1, 0.75 or 0.5 (flat, across the flow) |
| `gridDiametersUpstream` | 0 for no grid, otherwise the distance of the rods in bubble diameters |
| `particleDiameter_um`, `particleDensity_kg_m3` | particles sent past the bubble |
| `steps`, `warmup`, `seed` | time steps, steps before averaging starts, random seed |

Output: Re, Stokes number, number of averages, collisions on the front and rear halves, every
collision angle, and k and ε on a curve 1.5 times the bubble size around it. 20 000 steps take
about a minute and a half.

## Examples

`examples/collision-quartz.json`, `examples/rise-speed.json` and `examples/live-grid.json`.
The tutorials in `docs/tutorials/` walk through them.
