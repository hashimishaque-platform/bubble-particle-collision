# Contributing

Corrections, new models and new checks are welcome. Please open an issue first for anything
larger than a small fix, so the approach can be agreed.

## Ground rules

- **Every formula and every literature value has a source.** Give the full reference with
  DOI and the equation, table or page number, in the code comment and in the page's
  References.
- **Every new part has a check.** Add a row to `checks.js` that compares the new code with an
  exact solution or a published value. The check should fail (set `process.exitCode = 1`)
  when the result is out of tolerance.
- **No measured data without permission and a citation.** The program is a theoretical
  simulation. Measured data can only be added if it is published, its licence allows it and
  it is cited, and it then goes into a clearly marked comparison section.
- **Keep the style of the page:** plain HTML, Times New Roman, tables with grey borders,
  "Figure N." captions, short sentences, no frameworks and no web fonts. The page must work at
  400 px width without horizontal scrolling.

## How to work on it

```
node checks.js                  # all checks, about four minutes; rewrites checks.json
python3 build.py                # builds index.html from template.html and the .js files
node cli.js examples/rise-speed.json
```

The computation lives in plain JavaScript files without dependencies that run both in the
browser and in Node.js 18 or newer:

| File | Contents |
|---|---|
| `flow.js` | Axisymmetric Navier-Stokes solver, drag, particles, collision efficiency |
| `piv.js` | Synthetic PIV images and their evaluation |
| `swarm.js` | Bubble swarm: rise speed, shape, dipole flow, synthetic turbulence, tracks files |
| `lbm.js` | Live two-dimensional lattice Boltzmann simulation, particles, statistics |
| `template.html` | Text, figures and form of the page |
| `cli.js` | Batch runs without the browser |

Each module ends with a `module.exports` block for Node; `build.py` strips it when it puts the
page together. Do not edit `index.html` by hand.

Before a pull request: run `node checks.js` and `python3 build.py`, and commit the updated
`checks.json` and `index.html`. The workflow in `.github/workflows/checks.yml` repeats both
and fails if the committed files are not reproduced.
