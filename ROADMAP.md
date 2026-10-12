# Roadmap

Each step is developed on its own branch and merged into `main` once the checks pass
(`node checks.js`, which also rewrites `checks.json`). A step is ticked here with the date it was merged.
The project stays theoretical: no measured data and no data from any research centre.

## Responsiveness
- [ ] Run the flow solver in a Web Worker so the page stays responsive, with results identical to the current ones.
- [ ] Show the solver progress in Figure 1 while it runs.

## Generalised Sutherland equation
- [ ] Add Dukhin's generalised Sutherland equation for mobile surfaces, as given by Dai, Fornasiero and Ralston (2000), verified against numbers worked from the paper.
- [ ] Plot it with the inertial curve for clean bubbles and add a check row comparing the two.

## Partly mobile surface
- [ ] Mixed boundary condition for a stagnant cap: no slip for θ below the cap angle, shear-free above it; checks that cap angles of 0° and 180° reproduce the two existing cases.
- [ ] Drag against cap angle compared with published low-Re results (Sadhal and Johnson 1983), using values taken from the paper only.
- [ ] Form field for the cap angle and a figure of collision efficiency against cap angle.

## Attachment
- [ ] Sliding time of a particle along the bubble surface from the computed surface velocity.
- [ ] Attachment efficiency from an induction time (Yoon and Luttrell), and the collection efficiency as collision times attachment.

## Particle forces
- [ ] Basset history force as an option for heavy particles, with a check against an analytic settling case.

## PIV
- [ ] Window deformation as a third pass, and a comparison of the errors near the surface with the present method.
- [ ] Out-of-plane loss of particle pairs for a light sheet of finite thickness.
- [ ] Download the PIV vectors and the exact field as CSV files.

## Documentation
- [ ] `docs/numerics.md`: grid, deferred correction, boundary conditions and a grid-convergence table.
