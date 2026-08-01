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
pole lands at −1e-4 rad/s — a time constant of roughly three hours — so it
never tracks at all and the flow sits near zero. Bryson's rule with
honestly chosen maxima turns the same architecture into the best
controller in the set.

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

- **Standard PI** on flow error, Kp = 4000, Ki = 32000 (tuned on this
  plant; Ki/Kp = 8 places the PI zero on the mechanical pole at 8.2 rad/s).
  Anti-windup by conditional integration.
- **Feedforward PI**: adds α·V_ff with V_ff from inverse dynamics at the
  reference, assuming a *measured* load torque. α = 0.86 chosen from a
  sweep (rise time reaches its saturation-limited floor there; overshoot
  appears from α ≈ 0.90).
- **LQR + integral (Bryson)**: state augmented with ∫(Q_ref − y). Bryson
  maxima: 3 A current excursion, 40 rad/s speed error, 2.8e-5 m³/s·s
  integrated error (≈1 m³/h flow error corrected in 0.1 s), 24 V input.
  Closed-loop poles −40612, −78.9, −78.9.
- **LQR + integral (Q = I)**: identical architecture with the naive
  `Q = eye(3), R = 1` weighting, included to demonstrate the failure mode.
  Closed-loop poles −27302, −239, **−0.0001** — the integral state is
  effectively unweighted and the "tracker" degenerates to a regulator.

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

| Scenario | Controller | Rise [s] | Settle [s] | SS [m³/h] | Dip [%] | Recover [s] |
|---|---|---|---|---|---|---|
| Nominal | Standard PI | 0.111 | 0.308 | 85.00 | — | — |
| Nominal | Feedforward PI | 0.086 | 0.124 | 85.00 | — | — |
| Nominal | LQR (Bryson) | 0.086 | 0.115 | 85.00 | — | — |
| Nominal | LQR (Q = I) | never | >4 | 0.02 | — | — |
| Filter clogging | Standard PI | 0.111 | 0.308 | 85.02 | 0.76 | 0.113 |
| Filter clogging | Feedforward PI | 0.086 | 0.124 | 85.00 | 0.11 | 0 |
| Filter clogging | LQR (Bryson) | 0.086 | 0.115 | 85.00 | 0.26 | 0 |
| Filter clogging | LQR (Q = I) | never | >4 | 0.02 | n/a | n/a |

(Dip/recovery are undefined for a controller that never reached the
setpoint.)

Reading: feedforward PI and Bryson-weighted LQR both hit the
saturation-limited rise floor (~0.086 s) and never leave the ±0.5% band
under the load step; standard PI is ~30% slower to rise and takes 0.11 s to
recover. Feedforward PI shows the smallest dip, but that advantage rests
entirely on its assumption of a measured load torque — the LQR achieves
nearly the same with feedback only. The naive-weighted LQR fails outright,
which is the point of including it. No controller was tuned to make any of
them win; the PI baseline was tuned first, to be stable and
non-oscillatory, and left alone.

Figures: `figs/comparison.png` (all scenarios), `figs/zoom.png` (startup
transient and load-step dip, where the differences are actually visible),
`figs/ff_sweep.png` (feedforward scaling tradeoff).

## How to run

MATLAB (tested on R2025b) with the Control System Toolbox:

```
>> duct_fan_control
```

Regenerates the console output (`results.txt` content), both figures, and
the dt convergence check.

## Limitations

- LQR's classical margin guarantees do not survive input saturation; during
  the saturated startup the loop is effectively open.
- LQR requires full state feedback (current and speed measurements); PI
  needs only the flow output.
- The feedforward path assumes the load torque is measured; its disturbance
  rejection degrades with load-estimation error (see the α sweep).
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
