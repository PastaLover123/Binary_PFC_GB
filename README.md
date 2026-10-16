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

The supplied initialization/output code otherwise remains as provided,
including its duplicate bicrystal branch, smoothing, and movie-only output.
Composition must remain strictly between zero and one for the logarithms.
