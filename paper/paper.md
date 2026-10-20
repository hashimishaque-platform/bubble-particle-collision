---
title: 'Bubble and particle collision: an open, checked browser program for the flow around rising bubbles and particle collision in flotation'
tags:
  - JavaScript
  - flotation
  - bubble-particle collision
  - Navier-Stokes
  - lattice Boltzmann
  - particle image velocimetry
  - turbulence
authors:
  - name: Hashim Ishaque
    affiliation: 1
    # orcid: add the ORCID iD here before submission
affiliations:
  - name: add the affiliation here before submission
    index: 1
date: 20 October 2026
bibliography: paper.bib
---

# Summary

In froth flotation, mineral particles are collected by air bubbles rising through a
suspension. Before a particle can attach to a bubble it has to collide with it, and whether it
collides depends on the flow around the bubble and on the turbulence in the cell. *Bubble and
particle collision* is a single-page browser program, with a command-line interface, that
computes this chain of physics from published equations. It solves the steady, axisymmetric
Navier-Stokes equations around a bubble with an immobile or mobile surface, follows particles
with the Maxey-Riley equation [@Maxey1983], and computes the collision efficiency by
interception, settling and inertia, next to the correlations of Sutherland [@Sutherland1948],
Weber and Paddock [@Weber1983] and Yoon and Luttrell [@Yoon1989]. It makes synthetic PIV
recordings of the same flow and evaluates them as PIV software does, so every vector can be
compared with the exact flow. It draws bubble swarms with synthetic turbulence from a
kinematic simulation [@Fung1992], and it runs a live two-dimensional lattice Boltzmann
simulation [@Ginzburg2008] of a bubble in a turbulent stream that shows vortex shedding, the
turbulent kinetic energy and dissipation rate in the wake, and where particles hit the bubble.

# Statement of need

Collision models for flotation are usually used as closed-form correlations whose underlying
flow fields are not visible to the user, and comparisons between models, or between a model
and an experiment, require each group to reimplement the flow solution. Experimental work on
bubble-particle collision in turbulence, for example with particle tracking
[@Sommer2018], needs reference calculations with the same bubble size, Reynolds number and
particle Stokes number. This program provides them in a form that runs anywhere without
installation: the flow, the particle paths, the PIV evaluation and the turbulence statistics
are computed in the browser, every input can be shared as a link, and every figure can be
downloaded as data. The same code runs in Node.js for parameter studies.

The program is aimed at researchers and students in mineral processing and multiphase flow,
and at experimentalists who want to know how well PIV resolves the flow next to a bubble
surface before they build a set-up.

# Verification

Every part is checked against exact solutions or published values, and the checks run on
every change. The axisymmetric solver reproduces the drag, wake length and separation angle
of a sphere [@Johnson1999; @Clift1978] and the drag of a clean bubble [@Mei1994] for Reynolds
numbers from 0.1 to 150, on two grids. Collision efficiencies from integrated particle paths
agree with the exact value from the particle stream function to 0.07 %. The lattice Boltzmann
solver reproduces the decay of a Taylor-Green vortex to 0.1 % and the Strouhal number of the
cylinder benchmark at Re = 100 [@Schafer1996] (0.2985 against 0.295 to 0.305). The synthetic
turbulence has the kinetic energy of its spectrum to within its sampling error and no
divergence.

# Limitations

The axisymmetric solver is limited to steady flow (bubble diameters up to 1.2 mm with frother).
The bubble swarm is a kinematic picture without wakes. The lattice Boltzmann simulation is
two-dimensional, with an immobile surface. Attachment and detachment are not modelled.

# References
