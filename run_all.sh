#!/bin/sh
set -eu
cd "$(dirname "$0")"
python3 newton_demos.py > demo-output.txt
python3 - <<'PY' >> demo-output.txt
import json
d = json.load(open('values.json'))
g = d['gravity']
o = d['orbits']
print('Newton demos: all numerical tests passed')
print('10 s push: m1 %.3f m, m2 %.3f m (analytic 100, 50)' %
      (d['inertial_force']['masses']['1']['x_final'],
       d['inertial_force']['masses']['2']['x_final']))
print('Coasting: %.3f m at %.3f m/s; pulse final speed: %.3f m/s' %
      (d['inertial_force']['zero_force_coasting']['x_final'],
       d['inertial_force']['zero_force_coasting']['v_final'],
       d['inertial_force']['pulse_force']['v_final']))
print('Wrong force-as-speed endpoint: %.3f m; fixed pulse: %.3f m' %
      (d['inertial_force']['wrong_force_as_speed']['x_final'],
       d['inertial_force']['pulse_force']['x_final']))
print('Slope h=.01: %.5f; integral midpoint n=300: %.5f' %
      (d['calculus']['finite_difference_slopes']['0.01'],
       d['calculus']['integral_v_0_to_3']['300']['midpoint']))
print('Earth surface g: %.6f m/s^2' % g['earth_surface_acceleration_m_s2'])
print('Moon orbit: %.3f km/s, %.3f days' %
      (g['moon_orbit_speed_m_s'] / 1000, g['moon_orbit_period_days']))
print('Circular speed at 300 km: %.3f km/s' %
      (o['circular_speed_m_s'] / 1000))
for c in o['cases']:
    print('%g m/s -> %s, radius %.1f..%.1f km' %
          (c['speed_m_s'], c['outcome'], c['min_radius_m'] / 1000,
           c['max_radius_m'] / 1000))
print('Cold radial fall: %s' % o['zero_speed_radial_fall']['outcome'])
print('Euler circular energy drift: %.3g; Verlet: %.3g' %
      (o['coarse_euler_circular']['energy_drift_fraction'],
       o['corrected_verlet_circular']['energy_drift_fraction']))
for colour in ('red', 'violet'):
    ray = d['prism']['rays'][colour]
    print('Prism %s: %.3f deg deviation at 50 deg incidence' %
          (colour, ray['deviation_deg']))
print('Newton sqrt(2): %.15f' % d['newton_root']['root'])
PY
cat demo-output.txt
