# Numerical methods

This note describes how each part of the program is computed and how it is tested. The page
gives the equations; this note adds the grids, the boundary conditions and the convergence
data. All numbers in the tables come from `node checks.js`, which writes `checks.json`.

## 1. Flow around the bubble (Figures 1 to 4)

**Equations.** Steady, axisymmetric Navier-Stokes equations in the frame of the bubble, for
the Stokes stream function ψ and the vorticity variable ζ = ω r sin θ, made dimensionless
with the bubble radius a and the rise speed U (Re = U d / ν, with d the diameter).

**Grid.** Uniform in ξ = ln r and θ, so the spacing in r grows with the distance from the
bubble. The outer boundary is at R_max = 300 radii for Re < 1, 120 for Re < 10 and 60 above.

| Surface | Re | Grid "normal" | Grid "fine" |
|---|---|---|---|
| Immobile (no slip) | all | 90 points in θ, Δξ = 0.025 | 120 points in θ, Δξ = 0.018 |
| Mobile (shear-free) | below 20 | 90 points in θ, Δξ = 0.025 | 120 points in θ, Δξ = 0.018 |
| Mobile (shear-free) | 20 and above | 90 points in θ, Δξ = 0.0125 | 120 points in θ, Δξ = 0.009 |

A mobile surface has a thin boundary layer above Re of about 20, which needs the finer radial
spacing.

**Discretisation.** Second-order central differences. The convection terms use deferred
correction: first-order upwind differences in the matrix and the central difference as a
lagged source, so the converged solution is second order. The stream function is relaxed by
SOR (factor 1.8) and the vorticity by under-relaxation (factor 0.8, halved and restarted if
the iteration diverges). Iteration stops when the largest change of ζ per sweep, relative to
the largest |ζ|, is below 10⁻⁸ on the page and 10⁻⁹ in the checks.

**Boundary conditions.** On the surface ψ = 0, and either the velocity is zero (immobile
surface) or the shear stress is zero (mobile surface). At the outer boundary the flow is
uniform, ψ = r² sin²θ / 2 and ζ = 0. On the axis ψ = 0 and ζ = 0.

**Drag** follows from the surface vorticity and its normal derivative. **Wake length** and
**separation angle** are read from the sign change of ψ near the rear axis and of the surface
vorticity.

**Grid convergence.** The same cases on both grids, against published values:

| Case | Expected | Normal grid | Fine grid | Change, normal to fine |
|---|---|---|---|---|
| Immobile, drag coefficient, Re = 0.1 | 244.1 (Proudman and Pearson) | 244.7 | 244.5 | −0.08 % |
| Immobile, drag coefficient, Re = 50 | 1.57 (Johnson and Patel) | 1.567 | 1.571 | +0.26 % |
| Immobile, drag coefficient, Re = 100 | 1.08 | 1.068 | 1.077 | +0.84 % |
| Immobile, drag coefficient, Re = 150 | 0.90 | 0.860 | 0.871 | +1.3 % |
| Immobile, wake length / d, Re = 100 | 0.88 | 0.87 | 0.88 | |
| Immobile, separation angle, Re = 100 | 53.5° (Clift et al.) | 53.2° | 53.2° | 0 |
| Mobile, drag coefficient, Re = 10 | 2.477 (Mei et al.) | 2.442 | 2.441 | −0.04 % |
| Mobile, drag coefficient, Re = 100 | 0.375 | 0.377 | 0.376 | −0.27 % |
| Immobile, E / (d_p/d_b)², Re = 100 | 8.07 (Weber and Paddock) | 7.91 | 7.90 | −0.13 % |

The drag at Re = 150 changes most between the grids and is 4 % below Johnson and Patel on the
normal grid; use the fine grid near the upper end of the range.

## 2. Particles and collision efficiency

**Equation.** Maxey-Riley without history and lift forces, added-mass coefficient 1/2, Stokes
drag corrected with Schiller-Naumann at the settling speed (see the page).

**Integration.** Fourth-order Runge-Kutta, or below a Stokes number of 0.02 an exponential
integrator that treats the drag exactly (Cox and Matthews 2002). Velocities come from the
derivatives of the same interpolant of ψ that draws the streamlines.

**Collision efficiency.** Without inertia, exactly from the particle stream function
Ψ_p = ψ + (w_s/2) r² sin²θ at the grazing streamline. With inertia, by bisection on the
starting offset of paths that begin 8 radii upstream.

**Tests.** Exact limits for Stokes and potential flow; paths without inertia against the
exact value from Ψ_p (difference −0.07 % and −0.05 %); the two integrators against each other
for 30 µm quartz (difference −0.02 %).

## 3. Synthetic PIV

Tracers are moved half a time step back and forward along the computed flow. Particle images
are Gaussian spots integrated over the pixels, with background, noise and 8-bit quantisation.
Evaluation: FFT cross-correlation of mean-subtracted windows, three-point Gaussian sub-pixel
fit, a second pass with window offset, normalised median test (Westerweel and Scarano 2005).
Test: random error for a uniform shift, 0.037 px for 32 px windows.

## 4. Bubble swarm (Figures 5 and 6)

Rise speed from the Tomiyama et al. (1998) drag by bisection on the force balance; shape from
Wellek et al. (1966). The flow around each bubble is the potential-flow dipole of its
volume-equivalent sphere; the fields are added. The turbulence is a frozen kinematic
simulation (Fung et al. 1992): 200 Fourier modes with wavenumbers spaced geometrically from
2π/L to the smaller of 1/η and the finest scale of the figure grid, amplitudes
|A_n|² = |B_n|² = 2 E(k_n) Δk_n from the Kolmogorov spectrum with C = 1.62.

Tests: the surface condition u·n = U·n (error 3 × 10⁻⁹), zero divergence of both fields
(10⁻⁹ to 10⁻⁸ relative), the mean kinetic energy of the modes against the sum of E(k) Δk
(0.9991 over 40 000 points) and against the exact integral (1.000004), and the drag limits.

## 5. Live flow simulation (Figures 7 and 8)

**Method.** Lattice Boltzmann, D2Q9 lattice, two-relaxation-time collision (Ginzburg et al.
2008) with Λ = 1/4. The relaxation time follows from the lattice viscosity, τ = 3ν + 1/2.

**Grid and parameters.** 20 cells per bubble diameter; domain 300 × 120 cells (15 × 6
diameters); bubble centre 9 diameters from the inlet; lattice inflow speed 0.04 (Mach number
0.07). The viscosity follows from the Reynolds number, ν = 0.04 · 20 / Re.

**Subgrid model.** Smagorinsky with C = 0.15 in the form of Yu, Girimaji and Luo (2005):
τ = (τ₀ + (τ₀² + 36 C² |Q|)^½) / 2, with |Q| the norm of the non-equilibrium momentum flux.
It is always on in the page and in `cli.js`. The benchmark check below runs without it,
so it tests the underlying scheme.

**Boundary conditions.** Halfway bounce-back on the bubble and the rods (no slip); periodic
sides; velocity inlet after Zou and He (1997), with the inflow raised smoothly over the first
1500 steps; outlet of equilibrium populations at density 1 with the velocity of the cell
upstream, behind a damping layer of 40 cells where τ rises towards 1.

**Statistics.** Running averages every fourth step after 5000 steps: mean velocity,
k = (⟨u′²⟩ + ⟨v′²⟩)/2 and ε = 2ν(⟨s_ij s_ij⟩ − S_ij S_ij) with central differences.

**Tests.**

| Test | Expected | Result |
|---|---|---|
| Decaying Taylor-Green vortex, 64 × 64, 500 steps: amplitude / exp(−2νk²t) | 1 | 0.9989 |
| Cylinder in a channel at Re = 100 (Schäfer and Turek 1996, benchmark 2D-2), 20 cells per diameter: Strouhal number | 0.295 to 0.305 | 0.2985 |
| Particle released at rest in a uniform stream, after 3 response times: v / (U (1 − e⁻³)) | 1 | 1.000000 |

**Known limitations.** Two-dimensional; immobile surface; at Re of several hundred the
boundary layer is only a few cells thick, so the separation point and the smallest eddies are
approximate. The program caps Re at 1500.

## 6. Reproducing the numbers

```
node checks.js          # all checks, about four minutes; rewrites checks.json
python3 build.py        # rebuilds index.html from the sources and checks.json
```

The continuous-integration workflow in `.github/workflows/checks.yml` runs both on every push
and fails if `checks.json` or `index.html` differ from the committed files.
