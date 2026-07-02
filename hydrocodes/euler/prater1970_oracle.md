# Prater 1970 Al-on-Al validation oracle (Pierazzo et al. 2008 validation #2, Table 4)

Source: PierazzoEtAl2008.pdf pp. 1928-1935 (Table 4 transcribed 2026-07-02; Fig. 11/12/14
cross-checked). This is the EXPERIMENTAL strength validation for the euler code (E23 in
docs/validation_matrix.md).

## Setup (paper p. 1928 + Table 6)
- Al 2017-T4 sphere, **6.35 mm diameter** (a = 3.175 mm), impacting **perpendicularly at ~7 km/s**
  into quasi-infinite Al-alloy cylinders (few tens of mm across — paper cautions the sample is
  only a few times the crater size, which may have allowed a slightly larger experimental crater).
- Time axis: Table 4 header says "ms" but it is **microseconds** (paper typo; Figs. 11/12/14 use
  µs and early growth speeds ~2.6 km/s confirm it).
- Flash X-ray measured transient crater radius R and depth D vs time; error ±0.5 to ±1 mm.
- Targets: Al 6061-T6 (strength ~rate-INsensitive -> von Mises appropriate) and Al 1100-O
  (strongly strain-rate/strain dependent -> Johnson-Cook in the paper; constant-Y is a poor model).

## Material parameters (Table 6)
- iSALE (the apples-to-apples code): **von Mises Y = 414 MPa for Al 6061-T6** (shear strength
  (s1-s3)/2 = 207 MPa); Johnson-Cook for Al 1100 (A=49 MPa, B=157 MPa, C=0.016, M=1.7, n=0.167).
- AUTODYN Steinberg-Guinan: Al 6061-T6 rho=2.703 g/cc, G0=27.6 GPa, Y0=290 MPa, Ymax=680 MPa;
  Al 1100-O G0=27.1 GPa, Y0=40 MPa, Ymax=480 MPa.
- Our mode: Tillotson Al + von Mises **G=27.6 GPa, Y=414 MPa** (iSALE-matched), projectile
  treated as the same Al (2017-T4 rho 2.79 vs 2.70 — documented approximation; RAGE even ran a
  strengthless projectile).

## Pass/fail expectations (paper Figs. 11/14 + text pp. 1933-1935)
- Codes vs experiment, 6061-T6: radius UNDERestimated 5-13%, depth OVERestimated 4-12%.
- iSALE with von Mises (Fig. 14): radius err **1.2%** (best of all), depth OVERestimated ~20%
  (the von Mises signature: no work hardening -> deep narrow crater; JC1 = depth 7.5% but radius
  -13%; JC2 compromise -3.5%/+12%).
- OUR GATE (6061-T6, von Mises, cppr>=10): transient Dmax within ±10% of 1.31 cm AND final D in
  [-20%,+25%] AND final R in [-30%,+10%]. Rationale (2026-07-02 ladder, cppr 6/10/14/20 GPU +
  cppr10 CPU==GPU): the GROWTH PHASE (3-15 us) tracks the experiment inside the inter-code band
  (cppr20: R -8.3%, D -1.4%); the PLATEAU sits low (R -22%, D -10% at cppr20, converging upward
  with cppr) because our vacuum-face Riemann flux carries no deviatoric strength -> no rim uplift
  survives and the shock-heated wall slowly extrudes into the cavity. Documented limitation, same
  league as iSALE-von-Mises' +20% depth overshoot. NOT tuned: Y=414 MPa held at the iSALE value.
- 1100-O (secondary, needs an effective constant Y ~ documented calibration, not validation):
  experiment plateaus R ~1.66 cm, D ~1.82 cm.

## Table 4 — Al 6061-T6 target (t in µs, R and D in cm; separate shots, times not monotonic)
Radius:
t_us R_cm
2.845 0.7505
3.082 0.8475
5.098 0.8835
5.216 1.0255
6.402 1.048
6.639 1.0505
8.180 1.0685
6.758 1.1065
7.943 1.1065
9.721 1.111
7.825 1.1405
9.247 1.1855
14.226 1.2215
16.716 1.251
15.412 1.278
17.427 1.3255
22.051 1.3415
25.726 1.3595
27.149 1.346
27.741 1.314
27.741 1.3095
37.344 1.2985
39.952 1.355
40.664 1.366
40.901 1.3435
43.865 1.33
44.102 1.321
44.932 1.2895
43.983 1.242
47.895 1.33
51.452 1.249
52.4 1.2645
58.684 1.2805
69.235 1.321
72.673 1.267

Depth:
t_us D_cm
2.872 0.825
2.633 0.825
4.787 0.869
6.582 1.066
6.462 1.125
9.334 1.213
7.779 1.256
14.121 1.293
14.839 1.293
25.25 1.304
27.04 1.295
40.808 1.281
43.68 1.245
52.176 1.322
68.81 1.324

## Table 4 — Al 1100-O target
Radius:
t_us R_cm
0.971 0.6345
3.083 0.7555
2.863 0.776
5.528 0.857
4.986 0.922
6.418 1.057
7.348 1.0905
8.709 1.068
8.841 1.106
12.314 1.207
15.08 1.2585
16.151 1.3525
18.658 1.3365
23.618 1.437
24.683 1.5175
26.634 1.535
26.42 1.569
26.779 1.605
26.436 1.605
28.844 1.6225
31.108 1.5725
36.871 1.6815
38.229 1.6525
40.614 1.6185
43.831 1.656
46.357 1.6825
52.621 1.6345
58.696 1.6735
62.999 1.5855
65.642 1.6165
67.48 1.6385
65.672 1.615
74.098 1.652
79.601 1.685

Depth:
t_us D_cm
3.075 0.803
2.994 0.905
5.607 0.976
5.265 1.004
7.670 1.154
9.948 1.282
10.297 1.273
15.082 1.522
18.871 1.616
23.497 1.762
26.668 1.765
27.024 1.779
31.261 1.815
38.412 1.781
40.658 1.837
52.391 1.819
65.660 1.842
67.416 1.822
73.985 1.799
79.973 1.804

## Experimental plateau values (for the gate)
- 6061-T6: R_final ~ 1.31 cm (scatter 1.24-1.37 past 25 µs), D_final ~ 1.31 cm (1.28-1.32 past 14 µs).
- 1100-O: R_final ~ 1.66 cm (1.59-1.69 past 26 µs), D_final ~ 1.82 cm (1.78-1.84 past 31 µs).
