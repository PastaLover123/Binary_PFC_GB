# Binary PFC grain-boundary model

MATLAB GPU implementation of the supplied two-dimensional n–c XPFC model,
with Gaussian correlation kernels, inertial density and diffusive composition.

Run `simulateXPFC` after editing `genParams.m`. The simulation requires a
supported MATLAB GPU and Parallel Computing Toolbox; initialization uses
Image Processing Toolbox. The seed is triangular in 2D; kernel labels such
as BCC/FCC do not make the initializer three-dimensional.

## Dynamics

Density obeys
`n_tt + beta_n*n_t = M_n*laplacian(mu_n) + noise_n`.
Composition retains `c_t = M_c*laplacian(mu_c) + noise_c`.

The semi-implicit density update is

```text
n_k_next = ((2 + beta_n*dt)*n_k - n_k_prev + dt^2*NL_n)
           / (1 + beta_n*dt - dt^2*L_n)
L_n = -M_n*k^2
```

All division is pointwise. The driver initializes equal current and previous
density fields (zero discrete initial velocity). Save both `fields.n_k`
and `fields.n_k_prev`, together with composition and parameters, for a restart.
The existing driver does not implement checkpoint/restart loading.

The default damping is `beta_n = 1`; damping time is `1/beta_n` and long-time
density mobility is `M_n/beta_n`. The supplied `dt = 1` is retained as an
input setting, not a validated timestep for inertial dynamics. Check timestep
convergence for the chosen physical parameters. Free energy alone need not
decrease monotonically with density inertia.

The mixing contribution to the density chemical potential is positive;
accordingly `NL_n` contains `-F_mix` inside its positive `M_n*k^2` factor.

Gaussian kernels and the original interpolated correlation derivative are
retained. The latter assumes composition varies slowly on the density scale,
allowing it to pass through the density convolution derivative. No additional
composition filter or exact symmetrized derivative is introduced.

Noise amplitudes retain the original empirical normalization. Density noise
now acts as forcing in the second-order equation. Set both noise flags to
false for deterministic checks.

## Verification

Run `testDensityDynamics` in MATLAB for CPU-only checks of the actual step
routine: the linear density recurrence, density history, conservation of both
means, disabled noise, and the mixing-force sign.

The legacy initialization options retain their duplicate bicrystal branch
and smoothing. The seeded test below uses a separate initialization branch
and enables full-grid output.
Composition must remain strictly between zero and one for the logarithms.

## Seeded growth smoke test

Run `runDendriteTest` in MATLAB. This initializes a small triangular crystal
(radius 5 lattice units, smooth envelope width 1) in a uniform liquid on a
512 x 512 periodic grid. The density mean is exactly -0.05 and composition is
uniformly 0.1. Only the envelope is smoothed, preserving lattice oscillations.

The test disables both noises and movies, uses dt=0.01, beta_n=1, M_n=1,
M_c=0.01, and runs 20,000 steps (time 200). It displays density, composition
and strain every 1,000 steps and saves lossless full-grid MAT snapshots,
including density history, in a unique dendrite_test_* directory. These are
exploratory growth settings, not a calibrated dendritic-growth regime or a
validated stability/convergence result. A liquid seed in solid would test
melting; this test instead seeds solid in liquid to test solidification.

To continue longer, increase maxStep in runDendriteTest.m. Reduce dt if the
fields diverge or composition leaves (0,1); keep maxStep*dt fixed when
comparing timesteps. The driver stops before evaluating logarithms on invalid
composition. MATLAB GPU execution remains to be tested on the user's machine.

## Sigma5 shear test

In MATLAB:

```matlab
runSigma5ShearTest       % positive shear
runSigma5ShearTest(-1)   % reverse shear
runSigma5ShearTest(0)    % unstrained control, same duration
```

This uses a 2D **square** Sigma5 symmetric tilt bicrystal, orientations
+/-atan(1/3) = +/-18.434949 degrees, with grain 1 in the central half along Y
and grain 2 in the outer quarters. There are two periodic boundaries normal
to Y. The reference lengths are exactly 16*2*pi*sqrt(10)/k_ini in both
directions. A 405x405 odd grid gives about eight points per lattice unit
and avoids sheared Nyquist-plane ambiguity.

The test selects the existing two-Gaussian envelope at q1 and sqrt(2)*q1
for both species, under the explicit label "square". This changes the shell
selection from the single-peak dendrite test, while retaining the Gaussian
construction and the slowly varying composition approximation. Square
symmetry is seeded; its long-time stability must be checked in the control.
Initial mode amplitudes approximate the minimum within a two-shell square
ansatz, not a full binary equilibrium calculation. Composition begins uniform
at 0.1 and mean density at -0.05; no extra dislocations or noise are added.

The system relaxes for 5,000 steps (time 50), then applies simple shear
x=X+gamma*Y in increments every 100 steps. The signed rate is +/-5e-4 per
unit time, capped at |gamma|=0.1. dt=0.01, beta_n=1, M_n=1 and M_c=0.01.
After 25,000 total steps the loading reaches the cap. Both kernels and the
density/composition Laplacians use k^2=Kx^2+(Ky-gamma*Kx)^2. Fields and density
history stay on the material grid; shear changes the cell, not the array
values by an extra interpolation.

Each run has a unique sigma5_shear_* output directory. Full double-precision
fields and density history are saved every 1,000 steps. gb_tracking.csv and
gb_tracking.mat record time, strain, the two wrapped Y positions, unwrapped
displacements, signed velocities, and individual speed magnitudes. Velocities
are separate trailing linear fits over up to five output samples; no averaging
or velocity normalization is applied. Since y=Y for this shear, the measured
normal positions are also physical y positions. Live density images remain
in material coordinates; they do not draw the tilted cell outline.

The orientation tracker uses Fourier windows around each grain's first shell.
It marks tracking invalid rather than choosing arbitrary crossings if it no
longer resolves exactly one rising and one falling boundary. After a lost
identity, displacement/velocity values remain undefined for that run. Unwrapping
assumes each boundary travels less than half the box between output samples.

Independent Python checks passed for commensurability, initial crossings,
a known boundary translation, conjugate symmetry of sheared operators, and
200 nonlinear steps with conserved means and bounded composition. These
checks do not validate the full MATLAB GPU run or establish shear-driven
migration. Compare positive/negative loading with the zero-shear control and
check timestep and box-size convergence before interpreting velocities.

## 1024 x 1024 disconnection-pair test

```matlab
runSigma5DisconnectionTest             % c=0.3, pair, positive shear
runSigma5DisconnectionTest(0)          % pair, no shear
runSigma5DisconnectionTest(-1)         % pair, reverse shear
runSigma5DisconnectionTest(1,0.3,false) % same large box, flat GB control
runSigma5DisconnectionTest(1,0.1)       % lower composition, otherwise same test
```

The box contains 40 commensurate Sigma5 rectangular periods in each direction.
Its 1024 x 1024 grid resolves each reference lattice unit by about 8.10 points,
so this enlarges the physical box rather than only oversampling the small test.
The reference wave number is estimated with the actual interpolation weight
fractionB=3*c^2-2*c^3. Mean composition defaults to 0.3, and composition
mobility after release is 0.1 (previously 0.01). Increasing mobility changes
kinetics, not thermodynamic driving; omega and the Gaussian kernels are unchanged.

### Even-grid shear

Odd dimensions are not required. The operator builder now defines a mask that
removes only the Nyquist row and column on an even sheared grid, starting at
zero strain. Initialization projects both n and c into that subspace and
reconstructs their real fields. Each step projects the new density, its history,
and evolving composition into the same subspace. The zero Fourier mode is
retained. The ambiguous signed Nyquist vectors are therefore never evolved;
all retained modes use the original sheared radial operators. This is not a
general nonlinear dealiasing filter.

The driver reports the fraction of initial density spectral power discarded.
For this default seed it was about 1.22e-12 in the independent Python check.
That small value is specific to this resolved initialization, not a guarantee
for arbitrary fields. Simply discarding the imaginary part of an inconsistent
inverse transform would not provide the same consistent Fourier history.
The 405-square small test is still supported; 1024 is used for the large test.

### Seeded disconnection geometry

A finite terrace is placed on the lower boundary, between X=-Lx/8 and +Lx/8.
Its two ends have opposite step and Burgers content. In reference lattice
units a=2*pi/k_ini, the chosen square-Sigma5 mode is

```text
h = 3*a/sqrt(10)
b = 2*a/sqrt(10) along X
b/h = 2/3
```

This choice follows directly from the grain lattice translations
t_plus=(-a,3a)/sqrt(10) and t_minus=(a,3a)/sqrt(10).
Their difference is (b,0). Translating the boundary by t_minus and shifting
grain 1 by b reproduces the original two-grain registry:
n1(r+t_minus-b*eX)=n1(r), n2(r+t_minus)=n2(r).
This supplies a registry shift as well as a geometric step.

The terrace and phase shift use a periodic smooth window with core width
0.5*a. Grain 1's tangential displacement is localized about the lower GB,
decaying over 4*a into the bulk. The upper GB is initially flat. The density
mean is restored after constructing the seed. This is an approximate
disconnection-dipole seed, not a relaxed elastic solution or a Burgers-circuit
measurement of the relaxed cores. It can relax into another mode or annihilate.

### Relaxation and inspection

| Steps | Fields that evolve | Loading |
| --- | --- | --- |
| 1-5000 | Density only; c and c_k are held exactly fixed | None |
| 5001-10000 | Density and diffusive composition | None |
| 10001-30000 | Density and diffusive composition | Shear to +/-0.1 |

No noise is enabled. Density history is carried continuously through both
stage transitions. Initial and every-1000-step full-resolution MAT files are
saved; steps 5000 and 10000 therefore preserve the two relaxed states.

Figure 4 shows both local boundary contours Y(X), independently of the
spatially averaged GB trajectories. The contours are stored as status.localGB
and status.localGBx in each snapshot, and the latest contour is also saved in
gb_tracking.mat. Inspect these at steps 0, 5000 and 10000 to see whether the
terrace survives before interpreting a loaded run as motion of a seeded pair.
The tracker returns NaNs where crossings are ambiguous rather than inventing
a boundary. It does not automatically identify or count relaxed core types.
field_diagnostics.csv records both means, ranges, and composition standard
deviation; gb_tracking.csv keeps the two GB velocities separate.

For custom timing or parameters:

```matlab
p = makeSigma5DisconnectionParams(1,0.3,true);
p.freezeCompositionSteps = 2000;
p.shear.relaxSteps = 5000;
p.maxStep = 25000;
p.output_dir = fullfile(pwd,'my_sigma5_test'); % use a fresh folder
simulateXPFC(p);
```

Tests: testDensityDynamics exercises the actual CPU stepping routine, including
frozen composition and projected even-grid shear. testSigma5DisconnectionSetup
requires MATLAB GPU and checks the full initializer, lattice translations,
terrace height, mass conservation, and one frozen-composition sheared step.
MATLAB/Octave are unavailable in the development environment. The independent
Python translation passed initialization and 30 full-size nonlinear steps
(20 density-only, 10 coupled under shear), with conserved means, exact frozen
composition, bounded released composition and real inverse transforms.
This short check does not establish survival over the full relaxation stages
or migration of the relaxed pair.
