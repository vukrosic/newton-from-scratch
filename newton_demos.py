#!/usr/bin/env python3
"""Small, reproducible numerical demonstrations for a Newton film."""
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RUNS = ROOT / "runs"
G = 6.67430e-11
MU_EARTH = 3.986004418e14
R_EARTH = 6_371_000.0
R_MOON = 384_400_000.0

def inertial_force():
    dt, steps, force = 0.01, 1000, 2.0
    out = {"dt": dt, "steps": steps, "duration_s": steps * dt,
           "force": force, "masses": {}}
    for mass in (1.0, 2.0):
        x = v = 0.0
        samples = [{"t": 0.0, "x": x, "v": v}]
        for i in range(1, steps + 1):
            v += force / mass * dt
            x += v * dt
            if i % 100 == 0:
                samples.append({"t": i * dt, "x": x, "v": v})
        out["masses"][str(int(mass))] = {"samples": samples,
                                         "x_final": x, "v_final": v}
    x, v = 0.0, 2.0
    coasting = [{"t": 0.0, "x": x, "v": v}]
    for i in range(1, steps + 1):
        x += v * dt
        if i % 100 == 0:
            coasting.append({"t": i * dt, "x": x, "v": v})
    coasting_x = x
    x, v = 0.0, 0.0
    pulse = [{"t": 0.0, "x": x, "v": v}]
    for i in range(1, steps + 1):
        if i <= 200:
            v += force / 1.0 * dt
        x += v * dt
        if i % 100 == 0:
            pulse.append({"t": i * dt, "x": x, "v": v})
    out["zero_force_coasting"] = {"initial_v": 2.0, "samples": coasting,
                                   "x_final": coasting_x,
                                   "v_final": 2.0}
    out["pulse_force"] = {"force_duration_s": 2.0, "samples": pulse,
                          "x_final": x, "v_final": v}
    x, v = 0.0, 0.0
    wrong = [{"t": 0.0, "x": x, "v": v}]
    for i in range(1, steps + 1):
        if i <= 200:
            v = force / 1.0
        else:
            v = 0.0
        x += v * dt
        if i % 100 == 0:
            wrong.append({"t": i * dt, "x": x, "v": v})
    out["wrong_force_as_speed"] = {
        "label": "wrong illustrative model: force treated as speed",
        "force_duration_s": 2.0, "samples": wrong, "x_final": x,
        "v_final": v}
    duration = steps * dt
    out["analytic_x"] = {"1": force * duration**2 / 2,
                          "2": force * duration**2 / 4}
    return out

def impulse_pair():
    impulse, m1, m2 = 6.0, 2.0, 3.0
    v1 = v2 = 0.0
    before = {"a": m1 * v1, "b": m2 * v2}
    v1 += impulse / m1
    v2 -= impulse / m2
    after = {"a": m1 * v1, "b": m2 * v2}
    return {"impulse_on_a": impulse, "impulse_on_b": -impulse,
            "before_momentum": before, "after_momentum": after,
            "delta_momentum": {k: after[k] - before[k] for k in before},
            "total_momentum": sum(after.values())}

def calculus():
    slopes = {}
    for h in (1.0, 0.1, 0.01):
        slopes[str(h)] = ((3 + h) ** 2 - 3**2) / h
    integ = {}
    for n in (3, 30, 300):
        dt = 3.0 / n
        left = sum(2 * i * dt * dt for i in range(n))
        right = sum(2 * (i + 1) * dt * dt for i in range(n))
        mid = sum(2 * (i + 0.5) * dt * dt for i in range(n))
        integ[str(n)] = {"left": left, "right": right, "midpoint": mid,
                         "midpoint_error": mid - 9.0}
    return {"function": "x=t^2", "t": 3, "exact_slope": 6.0,
            "finite_difference_slopes": slopes, "integral_v_0_to_3": integ,
            "exact_displacement": 9.0}

def gravity():
    g_surface = MU_EARTH / R_EARTH**2
    g_moon = MU_EARTH / R_MOON**2
    speed = math.sqrt(MU_EARTH / R_MOON)
    period = 2 * math.pi * R_MOON / speed
    return {"G": G, "earth_radius_m": R_EARTH,
            "earth_mu_m3_s2": MU_EARTH,
            "moon_orbit_distance_m": R_MOON,
            "earth_surface_acceleration_m_s2": g_surface,
            "moon_orbit_acceleration_m_s2": g_moon,
            "moon_orbit_speed_m_s": speed,
            "moon_orbit_period_days": period / 86400}

def orbit(speed0, dt=2.0, verlet=True, sample_every=None):
    if sample_every is None:
        sample_every = 5 if dt == 2.0 else 1
    r0 = R_EARTH + 300_000.0
    period = 2 * math.pi * math.sqrt(r0**3 / MU_EARTH)
    limit = 2 * period
    x, y, vx, vy = r0, 0.0, 0.0, speed0
    t, min_r, max_r = 0.0, r0, r0
    samples = []
    energy0 = 0.5 * (vx * vx + vy * vy) - MU_EARTH / r0
    angular0 = x * vy - y * vx
    while t < limit and math.hypot(x, y) >= R_EARTH:
        if len(samples) == 0 or round(t / dt) % sample_every == 0:
            samples.append({"t": t, "x": x, "y": y,
                            "radius": math.hypot(x, y)})
        r = math.hypot(x, y)
        ax, ay = -MU_EARTH * x / r**3, -MU_EARTH * y / r**3
        if verlet:
            xn = x + vx * dt + 0.5 * ax * dt**2
            yn = y + vy * dt + 0.5 * ay * dt**2
            rn = math.hypot(xn, yn)
            axn = -MU_EARTH * xn / rn**3
            ayn = -MU_EARTH * yn / rn**3
            vx += 0.5 * (ax + axn) * dt
            vy += 0.5 * (ay + ayn) * dt
            x, y = xn, yn
        else:
            x += vx * dt
            y += vy * dt
            vx += ax * dt
            vy += ay * dt
        t += dt
        r = math.hypot(x, y)
        min_r, max_r = min(min_r, r), max(max_r, r)
    final_r = math.hypot(x, y)
    if not samples or samples[-1]["t"] != t:
        samples.append({"t": t, "x": x, "y": y, "radius": final_r})
    e = 0.5 * (vx * vx + vy * vy) - MU_EARTH / max(final_r, 1.0)
    angular = x * vy - y * vx
    circular_speed = math.sqrt(MU_EARTH / r0)
    if final_r < R_EARTH:
        outcome = "earth_contact"
    elif speed0 > math.sqrt(2.0) * circular_speed:
        outcome = "escape"
    elif abs(speed0 - circular_speed) < 1e-6:
        outcome = "circular"
    else:
        outcome = "elliptical"
    return {"speed_m_s": speed0,
            "integrator": "verlet" if verlet else "euler",
            "dt_s": dt, "samples": samples, "min_radius_m": min_r,
            "max_radius_m": max_r, "duration_s": t,
            "outcome": outcome,
            "observation_window": "two initial circular periods",
            "initial_specific_energy_j_kg": energy0,
            "final_specific_energy_j_kg": e,
            "energy_drift_fraction": (e - energy0) / abs(energy0),
            "angular_momentum_drift_fraction":
            (angular - angular0) / max(abs(angular0), 1.0)}

def orbits():
    r0 = R_EARTH + 300_000.0
    circular = math.sqrt(MU_EARTH / r0)
    cases = [orbit(s) for s in (3000.0, 6000.0, circular, 10000.0,
                                12000.0)]
    return {"radius_m": r0, "circular_speed_m_s": circular, "cases": cases,
            "zero_speed_radial_fall": orbit(0.0),
            "verlet_dt2": orbit(circular, dt=2.0, verlet=True),
            "coarse_euler_circular": orbit(circular, dt=20.0, verlet=False),
            "corrected_verlet_circular": orbit(circular, dt=20.0,
                                                verlet=True),
            "verlet_dt20": orbit(circular, dt=20.0, verlet=True)}

def prism():
    apex = math.radians(60.0)
    incidence = math.radians(50.0)
    rays = {}
    for colour, index in (("red", 1.51), ("violet", 1.53)):
        inside = math.asin(math.sin(incidence) / index)
        second_inside = apex - inside
        exit = math.asin(index * math.sin(second_inside))
        deviation = incidence + exit - apex
        minimum = 2 * math.asin(index * math.sin(apex / 2)) - apex
        rays[colour] = {"refractive_index": index,
                        "entry_angle_deg": 50.0,
                        "inside_angle_deg": math.degrees(inside),
                        "exit_angle_deg": math.degrees(exit),
                        "deviation_deg": math.degrees(deviation),
                        "minimum_deviation_bonus_deg": math.degrees(minimum)}
    return {"apex_angle_deg": 60.0, "common_incidence_deg": 50.0,
            "note": "illustrative chosen indices, not historical glass data",
            "rays": rays}

def newton_root():
    x, steps = 1.0, []
    for _ in range(8):
        fx, dfx = x * x - 2, 2 * x
        nxt = x - fx / dfx
        steps.append({"x": x, "f": fx, "next_x": nxt})
        x = nxt
    return {"target": "sqrt(2)", "steps": steps, "root": x,
            "absolute_error": abs(x - math.sqrt(2)),
            "derivative_zero_pitfall": {"x": 0.0, "f": -2.0,
            "derivative": 0.0, "result": "stop: tangent slope is zero"}}

def make_all():
    RUNS.mkdir(exist_ok=True)
    data = {"inertial_force": inertial_force(),
            "impulse_pair": impulse_pair(),
            "calculus": calculus(), "gravity": gravity(), "orbits": orbits(),
            "prism": prism(), "newton_root": newton_root()}
    for name, value in data.items():
        path = RUNS / (name + ".json")
        path.write_text(json.dumps(value, indent=2) + "\n")
    (ROOT / "values.json").write_text(json.dumps(data, indent=2) + "\n")
    return data

def tests(data):
    motion = data["inertial_force"]
    assert motion["duration_s"] == 10.0
    assert abs(motion["masses"]["1"]["x_final"] - 100.1) < 1e-9
    assert abs(motion["masses"]["2"]["x_final"] - 50.05) < 1e-9
    assert abs(motion["masses"]["1"]["v_final"] - 20.0) < 1e-12
    assert abs(motion["masses"]["2"]["v_final"] - 10.0) < 1e-12
    assert abs(motion["analytic_x"]["1"] - 100.0) < 1e-12
    assert abs(motion["analytic_x"]["2"] - 50.0) < 1e-12
    assert abs(motion["zero_force_coasting"]["x_final"] - 20.0) < 1e-12
    assert abs(motion["zero_force_coasting"]["v_final"] - 2.0) < 1e-12
    assert abs(motion["pulse_force"]["v_final"] - 4.0) < 1e-12
    assert abs(motion["wrong_force_as_speed"]["x_final"] - 4.0) < 1e-12
    assert abs(motion["wrong_force_as_speed"]["v_final"]) < 1e-12
    assert abs(data["impulse_pair"]["total_momentum"]) < 1e-12
    assert abs(data["calculus"]["finite_difference_slopes"]["0.01"] - 6) < .02
    estimate = data["calculus"]["integral_v_0_to_3"]["300"]["midpoint"]
    assert abs(estimate - 9) < 1e-10
    area = data["calculus"]["integral_v_0_to_3"]
    assert abs(area["300"]["left"] - 9) < abs(area["3"]["left"] - 9)
    assert abs(area["300"]["right"] - 9) < abs(area["3"]["right"] - 9)
    assert 9.7 < data["gravity"]["earth_surface_acceleration_m_s2"] < 9.9
    assert data["orbits"]["cases"][0]["outcome"] == "earth_contact"
    assert data["orbits"]["cases"][2]["outcome"] == "circular"
    assert data["orbits"]["zero_speed_radial_fall"]["outcome"] == \
        "earth_contact"
    circular = data["orbits"]["cases"][2]
    assert abs(circular["max_radius_m"] - data["orbits"]["radius_m"]) < 1e3
    assert data["orbits"]["cases"][4]["outcome"] == "escape"
    assert data["orbits"]["cases"][4]["final_specific_energy_j_kg"] > 0
    assert data["newton_root"]["absolute_error"] < 1e-12
    e = data["orbits"]["coarse_euler_circular"]["energy_drift_fraction"]
    v = data["orbits"]["corrected_verlet_circular"]["energy_drift_fraction"]
    assert abs(v) < abs(e)
    d2 = data["orbits"]["verlet_dt2"]["max_radius_m"]
    d20 = data["orbits"]["verlet_dt20"]["max_radius_m"]
    r0 = data["orbits"]["radius_m"]
    assert abs(d2 - r0) < abs(d20 - r0)
    assert abs(data["prism"]["rays"]["violet"]["deviation_deg"]) > \
        abs(data["prism"]["rays"]["red"]["deviation_deg"])

if __name__ == "__main__":
    values = make_all()
    tests(values)
    print(json.dumps({"status": "ok", "tests": "passed"}, indent=2))
