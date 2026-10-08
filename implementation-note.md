# Implementation note

`newton_demos.py` is the numerical source of truth. `values.json` is a
combined snapshot and `runs/*.json` are the same values split by scene so an
animation can load only the labels it needs. Running `run_all.sh` regenerates
all outputs deterministically.

The inertial example uses semi-implicit Euler integration for exactly 1000
steps and compares masses 1 and 2 against the analytic 10-second endpoints.
It also includes a 2 m/s zero-force coasting case and a 2 N pulse for 2 s,
then coasting. The impulse example shows equal and opposite momentum
changes. Calculus uses forward finite differences for `x=t^2` at `t=3`, plus
left, right, and midpoint rectangle estimates for `v=2t`.

The inertial JSON also contains a deliberately wrong comparison that assigns
`v = F/m` directly during the pulse. It is labelled as a dimensional mistake:
it stops at 4 m, while the integrated pulse reaches about 36.02 m.

Gravity uses `mu = 3.986004418e14 m^3/s^2`, Earth radius 6,371,000 m, Moon
mean distance 384,400,000 m, and `G = 6.67430e-11`. The Moon calculation is a
two-body point-mass orbit comparison.

Orbit cases use actual inverse-square acceleration. Velocity Verlet updates
position, then averages old and new acceleration for velocity. Integration
stops on Earth contact or after two initial circular periods. Samples are
every 5 steps at dt=2 s and every step at dt=20 s. Final states are always
included. A zero-speed radial fall, dt=2 versus dt=20 Verlet comparison, and
a coarse forward-Euler circular case are saved to expose numerical error.
Each orbit records specific energy and angular-momentum drift.

The prism uses one common 50-degree incidence angle, Snell refraction at both
faces, A=60 degrees, and chosen indices 1.51 (red) and 1.53 (violet). A
minimum-deviation calculation is retained as a separately labelled bonus;
the chosen indices are illustrative, not a claim about historical glass.
Newton's root demo applies tangent steps to `x^2 - 2` and records the
derivative-zero pitfall at x=0.
