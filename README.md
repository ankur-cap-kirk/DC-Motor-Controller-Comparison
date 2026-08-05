# Duct-fan airflow control: PI vs feedforward PI vs LQR

Controller comparison for a small ventilation duct fan driven by a brushed DC
micromotor. Three controllers regulate airflow to 85 m³/h at a 4000 rpm
operating point: a standard PI, the same PI with inverse-dynamics
feedforward, and an LQR with an augmented integral state weighted by
Bryson's rule. The study measures setpoint response, rejection of a
load-torque step (filter clogging), and behaviour under a continuously
varying load (gusty inlet).

The methodological point of the comparison: an integral-augmented LQR is
only competitive if the integral state is actually weighted. The three
states here live on wildly different scales (current ~1 A, speed ~400 rad/s,
integrated flow error ~1e-5 m³/s·s), so `Q = eye(3)` silently disables the
integral term. To make that concrete, the comparison includes the naive
`Q = eye(3), R = 1` design as a fourth trace: on this plant its integral
pole lands at −5.6e-5 rad/s — a time constant of about five hours — so it
never tracks at all and the flow sits near zero. Bryson's rule with
honestly chosen maxima turns the same architecture into the best
controller in the set.

Two further traces round out the set: an observer-based LQR that runs the
full-state design from the flow sensor alone (see *Does the LQR really
need two sensors?* below), and a pole-placement + static-prefilter design
with no integral state, kept as a second strawman (see *Why the integral
state?*).

## Hardware

**Motor — FAULHABER 2342 S 024 CR** (24 V, graphite commutation), all values
from the series 2342…CR datasheet, edition 2026-03-23:

| Parameter | Datasheet | SI value used |
|---|---|---|
| Terminal resistance R | 7.14 Ω | 7.14 Ω |
| Rotor inductance L | 264 µH | 2.64e-4 H |
| Torque constant Kt | 28.6 mNm/A | 0.0286 N·m/A |
| Speed constant kn | 334 rpm/V | Ke = 0.02859 V·s/rad |
| Rotor inertia | 5.8 g·cm² | 5.8e-7 kg·m² |
| Rated torque / current / speed | 20.5 mNm / 0.846 A / 5660 rpm | — |

Consistency check: in SI units Kt and Ke are the same physical constant;
here they agree to 0.03%, which validates the unit conversions.

**Load — 80 mm axial impeller, direct drive.** These are *estimates*, not
datasheet figures:

- Impeller inertia from mass and radius (20 g, 40 mm tip radius) using a
  solid-disc approximation J = ½mr². Conservative, since a real impeller
  carries much of its mass near the hub.
- Aerodynamic load assumed 8 mNm at the operating point (~3.4 W shaft
  power, plausible for an 80 mm fan at 4000 rpm), linearised as viscous
  drag about the operating point (secant), plus a viscous equivalent of the
  0.99 mNm datasheet friction torque.
- Airflow inferred from shaft speed by the fan affinity law (Q ∝ ω),
  calibrated at the nominal duct restriction: 85 m³/h at 4000 rpm.

Operating point vs ratings: 9.0 mNm, 0.314 A, 14.2 V — inside the 20.5 mNm
/ 0.846 A / 24 V continuous ratings with headroom for disturbances.

## Model and controllers

States x = [current; shaft speed], input u = terminal voltage (saturated to
[0, 24] V), disturbance τ_d = load-torque deviation from nominal [N·m],
output y = airflow [m³/s].

- **Standard PI** on flow error, Kp = 4002, Ki = 32076 — derived, not
  tuned: the electrical pole is 3300× faster than the mechanical one, so
  di/dt → 0 reduces the plant to first order (τ = 0.122 s,
  k_dc = 1.66e-3 m³/s per V, with the back-EMF damping Kt·Ke/R folded into
  the effective friction), and PI on a first-order plant gives an exactly
  second-order closed loop — two specs (ζ = 1.5, ωn = 20.9 rad/s, a
  deliberately overdamped baseline) determine both gains. The resulting
  Ki/Kp = 8 places the PI zero on the mechanical pole at 8.2 rad/s.
  Anti-windup by conditional integration.
- **Feedforward PI**: adds α·V_ff with V_ff from inverse dynamics at the
  reference, assuming a *measured* load torque. α = 0.86 chosen from a
  sweep (rise time reaches its saturation-limited floor there; overshoot
  appears from α ≈ 0.90).
- **LQR + integral (Bryson)**: state augmented with ∫(Q_ref − y). Bryson
  maxima: 3 A current excursion, 40 rad/s speed error, 2.8e-5 m³/s·s
  integrated error (≈1 m³/h flow error corrected in 0.1 s), 24 V input.
  Closed-loop poles −40612 and a well-damped conjugate pair
  −78.9 ± 39.4j (ζ ≈ 0.89, ωn ≈ 88 rad/s).
- **LQR + integral (Q = I)**: identical architecture with the naive
  `Q = eye(3), R = 1` weighting, included to demonstrate the failure mode.
  Closed-loop poles −27302, −239, **−5.6e-5** — the integral state is
  effectively unweighted and the "tracker" degenerates to a regulator.
- **LQR + observer**: the same Bryson LQR running on the flow measurement
  alone — a Luenberger observer (poles −2000, −20000) reconstructs the
  armature current. The observer is driven by the *delivered* (saturated)
  input and is given no knowledge of the load disturbance.
- **Pole placement + prefilter**: state feedback placed at the LQR's pole
  real parts (−40612, −78.9) with a static prefilter inverting the nominal
  DC gain — no integral state anywhere in the loop. It tracks exactly on
  the nominal plant and has no mechanism to notice the plant changed.

Simulation: forward Euler at dt = 2e-5 s. The stability guard uses the
fastest **closed-loop** pole (−40612 rad/s, LQR loop), not the open-loop
electrical pole; dt sits 2.4× under the Euler limit and halving it twice
moves the reported rise time by less than 0.2%. Disturbances fire at
t = 1.5 s, roughly 5× the slowest settling time, so setpoint metrics and
disturbance metrics cannot contaminate each other. `rng(0)` for
reproducibility.

## Results

From `results.txt` (rise = 10–90%, settle to ±2%, dip and recovery to
±0.5% for the +4 mNm filter-clogging step):

| Scenario | Controller | Rise [s] | OS [%] | Settle [s] | SS [m³/h] | Dip [%] | Recover [s] |
|---|---|---|---|---|---|---|---|
| Nominal | Standard PI | 0.111 | 0.00 | 0.308 | 85.00 | — | — |
| Nominal | Feedforward PI | 0.086 | 0.00 | 0.124 | 85.00 | — | — |
| Nominal | LQR (Bryson) | 0.086 | 0.04 | 0.115 | 85.00 | — | — |
| Nominal | LQR (Q = I) | never | n/a | >4 | 0.02 | — | — |
| Nominal | LQR + observer | 0.086 | 0.04 | 0.115 | 85.00 | — | — |
| Nominal | Pole place + prefilter | 0.085 | 0.00 | 0.114 | 85.00 | — | — |
| Filter clogging | Standard PI | 0.111 | 0.00 | 0.308 | 85.02 | 0.76 | 0.112 |
| Filter clogging | Feedforward PI | 0.086 | 0.00 | 0.124 | 85.00 | 0.11 | 0 |
| Filter clogging | LQR (Bryson) | 0.086 | 0.04 | 0.115 | 85.00 | 0.26 | 0 |
| Filter clogging | LQR (Q = I) | never | n/a | >4 | 0.02 | n/a | n/a |
| Filter clogging | LQR + observer | 0.086 | 0.04 | 0.115 | 85.00 | 0.02 | 0 |
| Filter clogging | Pole place + prefilter | 0.085 | 0.00 | 0.114 | 85.00 | 0.73 | never* |
| Gusty inlet | Standard PI | 0.113 | 0.43 | 0.315 | 85.06 | — | — |
| Gusty inlet | Feedforward PI | 0.087 | 0.06 | 0.126 | 85.01 | — | — |
| Gusty inlet | LQR (Bryson) | 0.087 | 0.15 | 0.116 | 85.00 | — | — |
| Gusty inlet | LQR (Q = I) | never | n/a | >4 | 1.17 | — | — |
| Gusty inlet | LQR + observer | 0.087 | 0.04 | 0.116 | 85.00 | — | — |
| Gusty inlet | Pole place + prefilter | 0.086 | 0.41 | 0.116 | 85.12 | — | — |

(Dip/recovery are undefined for a controller that never reached the
setpoint; overshoot likewise. \*The pole-placement design stays outside the
±0.5% band for the entire 1.25 s the load is applied — the offset is
permanent, not transient; see *Why the integral state?* below.)

Reading: feedforward PI and Bryson-weighted LQR both hit the
saturation-limited rise floor (~0.086 s) and never leave the ±0.5% band
under the load step; standard PI is ~30% slower to rise and takes 0.11 s to
recover. Under the continuously varying load the Bryson LQR holds the
tightest steady state (85.00 vs the PI's 85.06). The observer-based LQR
matches the full-state LQR everywhere while reading one sensor instead of
two. Pole placement + prefilter matches the LQR's rise on the nominal
plant, then parks 0.73% low for as long as the load persists. The
naive-weighted LQR fails outright, which is the point of including it. No
controller was tuned to make any of them win; the PI baseline was designed
first, to be stable and non-oscillatory, and left alone.

**How much of feedforward PI's edge is the measured-load assumption?**
Feedforward PI shows the smallest dip (0.11% vs the LQR's 0.26%), but it is
the only controller that gets told the load. Re-running the load step with
the feedforward seeing τ(1+ε) while the plant sees τ:

| ε | Dip [%] |
|---|---|
| −0.20 | 0.24 |
| 0.00 | 0.11 |
| +0.20 | 0.00 |

With the load *under*-measured by 20% the feedforward dip (0.24%) ties the
feedback-only LQR (0.26%) on depth — and the shapes favour the LQR: it
dips sharply but snaps back in ~0.1 s, while the under-measured
feedforward drags a slower recovery tail. Feedforward's margin over the
LQR is bought entirely with load knowledge; at 20% measurement error it
has none left:

![Feedforward with load-measurement error vs LQR](figs/ff_error.png)

### Does the LQR really need two sensors?

The full-state LQR reads armature current and shaft speed; the PI reads
only flow. Since (A,C) is observable — the rank check has been in the
script from the start — a Luenberger observer can reconstruct the current
from the flow measurement alone. With observer poles at −2000 and −20000
the estimated-state LQR matches full-state performance (0.086 s rise,
0.115 s settling) on a single sensor, so the extra-sensor cost was a
design choice, not a requirement.

Two implementation details carry the result: the observer is driven by the
input the actuator actually *delivered* (after saturation, not what the
controller demanded), and it has no E·τ term — it does not know the load
disturbance. That honesty has a curious consequence: under the load step
the observer-based LQR dips less than full-state (0.02% vs 0.26%). This is
not magic — the unmodeled load appears in the innovation, which steers the
current estimate the way a disturbance estimator would, and the controller
leans on the biased estimate in the helpful direction.

Observer pole placement is not free, but its cost only shows where the
estimate can actually be wrong. On the undisturbed startup the observer
has the exact model, the same initial state, and the same delivered input,
so the innovation is identically zero and *any* pole choice reproduces
full-state behaviour — which is why the sweep below runs the
filter-clogging step, the one scenario with an estimation-error source:

| Observer poles | Dip [%] | Recover to ±0.5% [s] |
|---|---|---|
| −300, −3000 | 0.05 | 0.129 |
| −800, −8000 | 0.01 | 0.057 |
| **−2000, −20000** | **0.02** | **0** |
| −300, −60918 | 4.42 | 0.053 |

Too slow and the estimate lags the disturbance, stretching recovery past
the standard PI's 0.112 s; badly split — one pole slower than the
controller's dominant pair, one 1.5× its fastest — and the slow estimator
mode turns a sub-0.3% event into a 4.4% excursion. The chosen pair is the
only one that never leaves the band. There is a window, and locating it is
part of the design.

### Why the integral state? (pole placement + prefilter)

Pole placement at the same dominant speeds with a static prefilter
(u = −K·x + G·Q_ref, G inverting the nominal DC gain) tracks perfectly on
the nominal plant: 85.000 m³/h, 0.085 s rise, no integral state needed.
Then the filter clogs, and it settles at **84.378 m³/h — a permanent 0.73%
error** for as long as the load persists. A static prefilter inverts the
nominal plant once, at design time, so it has no mechanism to notice the
plant changed; the integral state is that mechanism. This is the same
lesson the naive `Q = eye(3)` trace teaches about weighting, one level up
the architecture: both traces are correct designs by their own textbook
recipe, and both fail for a reason the winning design explicitly handles.

### Startup and load-step response (zoom)

The naive LQR runs along the floor of the left panel; the right panel is
where the tracking controllers actually differ — note the pole-placement
trace holding its permanent offset at 84.38 while everything with an
integral state climbs back.

![Startup transient and load-step dip](figs/zoom.png)

### All scenarios

![Controller comparison, all scenarios](figs/comparison.png)

### Feedforward scaling tradeoff

Rise time reaches its saturation-limited floor by α ≈ 0.86; overshoot
appears from α ≈ 0.90. α = 0.86 is the knee.

![Feedforward scaling sweep](figs/ff_sweep.png)

## How to run

MATLAB (tested on R2025b) with the Control System Toolbox:

```
>> duct_fan_control
```

Regenerates the console output (`results.txt` content), all four figures,
and the dt convergence check.

## Limitations

- LQR's classical margin guarantees do not survive input saturation; during
  the saturated startup the loop is effectively open.
- The PI baseline was designed conservatively (ζ = 1.5) for an
  unambiguously non-oscillatory response, then held fixed. The gain grid in
  the script shows a more aggressive pairing at the same Ki/Kp ratio
  (Kp = 8000, Ki = 64000) reaches the same 0.086 s rise time as the LQR, so
  the reported rise-time gap reflects this choice rather than a structural
  advantage. The LQR's advantages in settling time, disturbance dip, and
  behaviour under naive weighting are unaffected.
- The feedforward path assumes the load torque is measured. The α sweep
  tests feedforward *authority*, not measurement error; the τ(1+ε) test in
  Results tests the measurement assumption directly, and ±20% error erases
  feedforward's advantage over the LQR.
- The mechanical dynamics rest on two estimates: impeller inertia is 28×
  the rotor inertia, and the damping comes from the assumed aero load. The
  sensitivity study in the script re-runs the load-step scenario with J and
  Bd each at ±50% (controllers fixed at their nominal designs). The numbers
  move — rise times span 0.041–0.161 s across the four corners — but the
  ranking never changes: feedforward PI and LQR stay ahead of standard PI,
  and the dip ordering FF (≤0.12%) < LQR (≤0.29%) < PI (≤0.82%) holds in
  every corner.
- Forward Euler integration, with the convergence check stated above.
- The speed→flow map is the affinity law at fixed nominal restriction;
  restriction changes are modelled only as shaft-torque disturbances, so
  "flow" during a disturbance is speed-referenced, not duct-measured.
- Impeller inertia and aerodynamic load are stated estimates, not
  measurements.

## References

Modelling, controller structure, and the weighting methodology follow
standard texts: Nise, *Control Systems Engineering*; Ogata, *Modern
Control Engineering*; Williams & Lawrence, *Linear State-Space Control
Systems*.
