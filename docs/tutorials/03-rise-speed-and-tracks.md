# Tutorial 3: rise speed, shape and a tracks file

**Question.** How fast does a bubble rise, how flat is it, and how does that change with a
surfactant? How do bubble tracks from an experiment compare?

## From the command line

```
node cli.js examples/rise-speed.json --csv rise.csv
```

For each diameter and surface tension (72.75, 50 and 30 mN/m) the file gives the Eötvös
number, the rise speed after Tomiyama et al. (1998) for pure, slightly and fully
contaminated water, and the aspect ratio after Wellek et al. (1966). A lower surface tension
raises Eo and lowers the rise speed of bubbles above about 2 mm.

## In the page

1. In "Bubble swarm in a flotation cell", choose "chain above one needle", 2.5 mm, and set the
   surface tension. Press Play: the bubbles rise and their flow moves with them.
2. Figure 6 shows the bubbles of Figure 5 as a red cross on the correlations.

## Your own tracks

Choose "tracks file" and load a CSV with the columns `TimeStep, ID, X, Y, Z, VX, VY, VZ, a, c,
alpha` (positions and semi-axes in mm, y upwards, velocities in m/s, alpha in degrees). The
file is read in your browser only. The page:

- estimates the time between steps from the distance moved and the velocity, and shows it
  as frames per second; check this against the recording;
- plays the time steps with the flow around each bubble;
- puts every track in Figure 6 at its mean size, rise speed and aspect ratio. Points above the
  curves in (a) and below the curve in (b) mean bubbles that behave as in clean water.

"Figure 6 data (CSV)" saves the curves and the track statistics together.
