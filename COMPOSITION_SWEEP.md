# Plain Sigma5 composition sweep

Pull main, then run in MATLAB:

```matlab
folder = runSigma5CompositionSweep;
```

The five default mean compositions are [0.1 0.2 0.3 0.4 0.5]. Cases run
sequentially on the GPU. No disconnections or cut-and-shift preparation are used.
No disconnection initializer is called.

The grid is 256 points along the GBs (X) and 1024 normal to them (Y), with 10
and 40 commensurate rectangular Sigma5 repeats respectively. There are about
8.1 points per reference lattice unit. The two boundaries begin half a Y box
apart, approximately 63.25 reference lattice units. This is a longer
separation and a shorter parallel length than the earlier small square test.

All compositions share the exact same physical box, reference density,
reference wave number (chosen at c_ref=0.3), initial lattice amplitudes, timestep
and imposed shear history. Species lattice constants and Gaussian kernels are
unchanged. Mean density is -0.05; M_n=1, beta_n=1, M_c=0.1 and dt=0.01.
Noise is off. Both density and composition evolve during the 2,000-step
unstrained relaxation. Simple shear then ramps at 5e-4 to strain 0.1 over
20,000 further steps. Total time is 220; loading lasts 200 time units.

Keeping the physical box fixed means composition-dependent lattice mismatch
can produce different residual stresses. This is a comparison under a common
loading protocol, not a separation of pure solute drag from all other
composition-dependent effects or a constant-stress mobility measurement.

## Outputs

One four-panel figure compares all five cases:
- signed displacement of GB 1;
- signed displacement of GB 2;
- speed magnitude of GB 1;
- speed magnitude of GB 2.

The baseline and time zero are the end of the unstrained relaxation. Velocities
are separate trailing fits over up to seven output samples (six time units at
the default cadence). Fits never span the initial relaxation/loading transition.
No averaging of the boundaries and no velocity normalization are applied.

The result directory contains composition_comparison.fig and a vector
composition_comparison.pdf, plus velocity_summary.csv. The summary velocities
are signed straight-line fits over strain 0.05-0.1, with separate magnitudes;
they are late-loading summaries, not claimed constant velocities.

Each case saves gb_tracking.csv, tracking.mat with its Y profiles, and initial,
end-of-relaxation and final full-resolution snapshots. The compact root
composition_sweep.mat checkpoint updates at every output. Snapshot saving can
be disabled. Only a 1D orientation profile is gathered during tracking; no
live whole-domain graphics or per-X contour detection slows the sweep.

To replot completed or interrupted results without rerunning:

```matlab
plotSigma5CompositionSweep(folder)
```

Or call it without arguments and select the result directory.

An unresolved orientation profile or grain width below six lattice units stops
that case; the remaining compositions still run. Failures and partial cases
are recorded in velocity_summary.csv. Missing identities are not interpolated
across gaps. Inspect the saved density/profile if a case is unresolved or
appears to transform away from the square bicrystal.

## Custom runs

```matlab
folder = runSigma5CompositionSweep([0.05 0.15 0.25 0.35 0.45]);
```

Optional short smoke test:

```matlab
opts = struct('relaxSteps',100,'loadingSteps',200,'saveSnapshots',false);
folder = runSigma5CompositionSweep([0.1 0.2 0.3 0.4 0.5],opts);
```

Durations must be multiples of 100. The short test does not reach the 0.05-0.1
summary interval, so its summary speeds will be undefined. For production
changes to grid, timestep, reference geometry or mobility, edit
makeSigma5CompositionParams.m and recheck convergence.

The simulations require the same supported MATLAB GPU as the previous tests.
testSigma5MotionAnalysis is a CPU-only test of separate velocities, periodic
unwrapping, loading baselines and missing-data handling. MATLAB is unavailable
in the development environment; numerical translation checks do not replace
running the MATLAB implementation on the user's GPU.
