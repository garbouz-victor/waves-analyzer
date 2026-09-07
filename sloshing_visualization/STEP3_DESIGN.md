# STEP 3A — design and pre-implementation self-review

Base: `c6dcf3faf28ff2cf9f303e5d5d8832d9ed2f7517`. This is a **new physical model**;
the linear solver, historical reports and existing FEM datasets are not inputs to
its equations. Only benchmarks are authorized here. No tank release or film movie
is authorized by this design's readiness status.

## 1. Dimensional model

Use the volume-averaged, divergence-free velocity of the
[Abels–Garcke–Grün formulation](https://arxiv.org/abs/1104.1336), with an explicitly
derived conservative-potential extension below. Coordinates are `(x,z)`, z upwards.
`phi=+1` denotes liquid; `phi=-1` denotes gas. Define

```
rho = rho_bar + rho_prime phi
rho_bar = (rho_l+rho_g)/2; rho_prime = (rho_l-rho_g)/2
f = lambda/(4 epsilon) (phi²-1)² + lambda epsilon/2 |grad phi|²
Psi = g z - a_x x;  b = -grad Psi = (a_x,-g)
mu_0 = lambda/epsilon (phi³-phi) - lambda epsilon Delta phi
mu_ch = mu_0 + rho_prime Psi
J = -rho_prime M grad(mu_ch)

div u = 0
phi_t + div(phi u) = div(M grad mu_ch)
(rho u)_t + div(u tensor (rho u+J))
    = -grad p + div(2 eta D(u)) + mu_0 grad phi + rho b
```

`p` is mechanical pressure (Pa), `eta` dynamic viscosity (Pa s), `M>0`
mobility (m²/(Pa s)), and `mu_ch` has units Pa. In index notation the extra
momentum transport is `partial_j(u_i J_j)`, not `partial_j(J_i u_j)`.
No assumption of small velocity or small displacement is made.

For the weak implementation use the exactly equivalent pressure
`pi=p-phi mu_0+rho_bar Psi`. The entire nonviscous force then becomes
`-grad pi-phi grad mu_ch`; **do not add rho b again**. Store/recover mechanical
pressure as `p=pi+phi(mu_ch-rho_prime Psi)-rho_bar Psi`. A Laplace test of `pi`
alone would be wrong. The gravity contribution to chemical potential is explicit:
omitting it while claiming decay of gravitational + interfacial energy leaves
the uncontrolled power `integral J dot grad(Psi)`.

The first benchmark set uses matched density or explicitly modest density ratios.
It is not water–air. Affine density is essential for the following mass identity;
silently clipping it would change that identity. Use affine viscosity initially,
check positivity of **both** coefficients including phase overshoot, and reject
inadmissible states. Large-density-ratio production is not qualified by this choice.

## 2. Energy normalization and law

For `phi(s)=tanh(s/(sqrt(2) epsilon))`, the two interfacial energy terms are
equal and direct integration gives `sigma = 2 sqrt(2) lambda/3`. Therefore
`lambda=3 sigma/(2 sqrt(2))`, not sigma itself. This will have independent
symbolic, quadrature and pressure-jump checks.

Per unit out-of-plane length (energies J/m, powers W/m),

```
E = integral_Omega [rho |u|²/2 + rho Psi + f] dA
      + integral_wetting_walls f_wall(phi) ds
dE/dt = - integral_Omega [2 eta D:D + M |grad mu_ch|²] dA
        - integral_slip_walls beta |u_tau|² ds
```

This is a continuum law for stationary impermeable walls, zero chemical flux,
fixed Psi and no inlet/outlet work. A release changes only `a_x`; the pre-release
and post-release potential-energy references must not be mistaken for dissipation.
Driven falling-film benchmarks instead require an external-work budget.

Derivation: the CH equation implies
`rho_t+div(rho u+J)=0`. This cancels the kinetic transport terms, including J.
Multiplying CH by mu_ch cancels the capillary work in momentum. The pressure does
no work for incompressible impermeable flow. Viscous stress and Navier friction
give the two remaining positive dissipations. Wall variation vanishes by the
variational boundary condition below. A conservative CH weak form preserves its
constant test function even with weakly, not pointwise, solenoidal Taylor–Hood u.

## 3. Walls and moving contact lines

Use the **instantaneous variational wetting law** in this first implementation:

```
f_wall = -sigma cos(theta_e) (3 phi-phi³)/4
lambda epsilon partial_n phi + f_wall'(phi) = 0
n dot M grad mu_ch = 0
u dot n = 0
beta u_tau + (2 eta D(u) n)_tau = 0; beta=eta/L_s
```

Here theta is measured **through the liquid**. The wall energies obey
`f_wall(-1)-f_wall(+1)=sigma cos(theta_e)`, Young's convention. The microscopic
equilibrium angle is fixed, but the apparent angle fitted several epsilon away
from the wall need not equal it during motion. There is no hysteresis.

This is not a finite-rate wall-relaxation model. In the notation of
[Qian–Wang–Sheng](https://arxiv.org/abs/cond-mat/0602293), the uncompensated wall
variation L is zero. Introducing finite relaxation would require both its
dissipation and the corresponding generalized Navier/Young stress; it must not
be bolted onto pure Navier slip alone. Wall-relaxation dissipation is identically
zero for the selected law, not an omitted positive contribution.

Diffuse interface + positive microscopic L_s permit contact motion without a
sharp/no-slip singularity. Neither parameter is automatically molecularly accurate.
Side walls are wetting/slipping; tank bottom would be no-slip, with its wetting
law explicit; tank lid would be impermeable/free-slip and nonwetting. Benchmark
boundary sets are stated separately, particularly the horizontal wetting wall
used to measure droplet angles. No wall film is prescribed or painted in.

## 4. Discretization and conservation review

Triangles, velocity CG2 vector / pressure CG1; phase and chemical potential use
the same CG degree, initially CG2 with a CG1 comparison. Use full symmetric D,
consistent quadrature, UFL automatic Jacobian and monolithic PETSc SNES. Closed
containers require a pressure gauge; it must remove only the nullspace, not add
compressibility. Report strong divergence separately from the pressure constraint.

Use backward-Euler startup then constant-step BDF2. The momentum transport is
written in skew form with mass flux `rho u+J`, and temporal counterpart
`D_t(rho u)-0.5 D_t(rho) u`. In the continuum this is the conservative AGG
equation after using mass balance. This helps avoid artificial transport power
from weak divergence. CH advection remains conservative. **BDF2 with a quartic
potential is not claimed unconditionally energy stable or to obey an exact
physical-energy identity.** Measure budget defects and refine dt. Failed Newton
steps are explicit; fixed-dt validation cannot silently change its step size.

### Self-review before implementation

- Solenoidal u is volume-averaged, not barycentric: compatible with unequal rho.
- Density is affine in phi, so the CH-to-density mass identity is exact.
- J is included in momentum in the correct tensor orientation.
- Conservative CH preserves integral phi; no fictitious mass correction is used.
- Gravity enters total chemical potential and the pressure transform consistently;
  it is not double-counted as a second momentum body force.
- Pressure reported for Laplace verification is p, not pi.
- Wetting variation and Navier work close the wall-energy budget for L=0.
- Positive rho/eta are checked, not assumed outside phi in [-1,1].
- Continuum consistency does not prove spatial, temporal or nonlinear convergence.

**No hidden algebraic conflict found in this continuum specification.** This is
not a numerical validation verdict; every implementation claim still needs tests.

## 5. Scales, mesh and computational budget

Each config specifies L and U (not inferred from an arbitrary transient maximum).
Report Re=rho_l U L/eta_l, Fr=U/sqrt(g L), We=rho_l U² L/sigma,
Bo=(rho_l-rho_g)g L²/sigma, Ca=eta_l U/sigma, Cn=epsilon/L,
S=L_s/L, density/viscosity ratios, and Pe_CH=U L²/(M sigma), using
mu-scale sigma/L. A zero gravity or density contrast makes the associated
capillary length infinite, not a divide-by-zero or a finite invented value.

`l_c=sqrt(sigma/(Delta rho g))`. Report a/l_c and epsilon/l_c. Define the diffuse
transition explicitly as phi=-0.9 to +0.9, whose planar width is
`2 sqrt(2) atanh(0.9) epsilon = 4.164... epsilon`. Require at least eight cells
across that width for qualification; also report epsilon/h_normal itself.
Refine independently in h, dt, epsilon, M and L_s. Grading must resolve wall-normal
films; report anisotropy and minimum angle, not h_min alone.

The local Docker daemon has 8 CPUs and 8.06 GB RAM; host free disk was 34 GB at
audit. A uniform n by n triangular grid with CG2 phase has approximately
`4(2n+1)²+(n+1)²` scalar unknowns: 17,989 at n=32; 70,789 at n=64;
280,837 at n=128. Sparse LU fill can dominate memory. Start well below 100k
unknowns, measure matrix nnz/memory and wall time before larger refinements.
No depth-10 tank mesh is approved on estimates alone.

## 6. Environment and ordered benchmark gates

Use official `dolfinx/dolfinx:v0.10.0`, Linux amd64, real double precision.
The immutable image digest and exact Python/PETSc/petsc4py/mpi4py/UFL/Basix/
NumPy/SciPy/h5py/Matplotlib versions will be recorded from the actual pulled image
in `docker/step3/stack.json` before solver execution. Use the 0.10
[SNES NonlinearProblem API](https://docs.fenicsproject.org/dolfinx/v0.10.0/python/generated/dolfinx.fem.petsc.html),
direct LU/MUMPS for small validation problems. No changes to the old Python venv.

Gates, in order:

1. Symbolic normalization, mass/energy/sign identities; mesh/parameter guards.
2. Flat stationary interface: mass, energy, spurious velocity vs h and phase order.
3. Circular droplets: mechanical pressure jump vs sigma/R, several R, h, epsilon.
4. Static 60/90/120 degree droplets: contour fits away from wall, fit-window study.
5. Spreading: correct movement, mass, energy, dt, slip and mobility sensitivity.
6. Falling film: velocity/flux against a declared analytic boundary-value problem.
   Navier slip adds `g L_s h²/nu` to the single-liquid no-slip flux; finite gas
   shear must be included in the reference rather than blamed on the solver.
7. Separate discretization and regularization studies; MPI n=1/n=2 check;
   restart reproduces BDF history and energy accounting.

Only measured benchmark results can authorize a bridge run. Lack of resources,
unfinished studies or a fundamental failure produce **MODEL NOT YET VALIDATED**,
not a relaxed threshold. Film thickness h<2 epsilon or <4 wall-normal cells is
unresolved; h>=4 epsilon and >=6 cells is only a resolution candidate, still
requiring model-sensitivity studies. No film claim follows from a phi heatmap.
