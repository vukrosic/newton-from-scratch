# Become Isaac Newton — from scratch

Small runnable experiments behind the animated lesson. Requires Python 3, with no third-party packages or network calls.

```sh
./run_all.sh
```

The launcher works from any directory. It runs numerical checks and regenerates runs/*.json, values.json and demo-output.txt.

Experiments: force and coasting; shrinking finite differences and accumulated area; inverse-square gravity; five launch outcomes; Euler drift versus Verlet; a separately integrated Moon orbit; illustrative prism refraction; Newton root-finding.

Every orbit is numerical solver output. The Earth is fixed, gravity is a point-mass field outside a spherical boundary, and there is no atmosphere. The Moon run is an ideal circular Earth-only comparison, not a precision ephemeris. Body sizes are enlarged in the film for readability. Prism indices are chosen illustrative values, not historical measurements.

Change launch speed, timestep or force duration in newton_demos.py and rerun. Compare timestep convergence, energy drift and observation windows before trusting a path. Real orbits add other bodies, drag, uncertainty and relativistic effects.

This repository contains the executable numerical demonstrations shown in the film. Animation and narration remain in the production project.
