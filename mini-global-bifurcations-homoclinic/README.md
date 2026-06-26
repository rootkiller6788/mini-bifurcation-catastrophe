# Global Bifurcations & Homoclinic Orbits

## Large-scale reorganizations of phase space structure

This module implements detection and analysis of global bifurcations —
dynamical events that cannot be captured by local linearization and involve
large-scale reorganizations of the phase portrait.

### Theoretical Foundation

While local bifurcations (saddle-node, Hopf, pitchfork) are detected by
eigenvalue crossings, **global bifurcations** involve the collision of
invariant manifolds and require nonlocal analysis:

- **Homoclinic bifurcation**: stable and unstable manifolds of a saddle collide
- **Heteroclinic cycle**: sequence of saddle connections forming a closed loop
- **Blue-sky catastrophe**: periodic orbit disappears without period blowup
- **Homoclinic tangency**: onset of chaotic dynamics via Smale horseshoe

### Shil'nikov Theorem (1965)

A homoclinic orbit to a saddle-focus with saddle quantity sigma > 0 implies
infinitely many periodic orbits and chaotic dynamics:

```
sigma = -Re(lambda_s) / Re(lambda_u)
```

If sigma > 1: stable periodic orbit. If sigma < 1: chaotic dynamics.

### Melnikov Method

For near-integrable perturbed systems, the Melnikov integral measures
the splitting of stable and unstable manifolds:

```
M(t0) = integral_{-inf}^{inf} f(q0(t)) ^ g(q0(t), t+t0) dt
```

Simple zeros of M(t0) indicate transverse homoclinic intersections
and the onset of chaotic dynamics via Smale's horseshoe construction.

### Key Methods

1. **Shooting method**: Integrate from unstable manifold, check return proximity
2. **Lin's method**: Boundary-value problem formulation for homoclinic orbits
3. **Melnikov method**: Perturbation-theoretic detection for near-integrable systems
4. **Poincare map analysis**: Reduced return map near homoclinic connection
5. **Continuation**: Track homoclinic orbits in parameter space

### Implementation

| File | Purpose |
|------|---------|
| `homoclinic_core.h/c` | Shooting-based homoclinic orbit detection |
| `global_bifurcation.h/c` | Bifurcation detection (homoclinic, heteroclinic, blue-sky) |
| `melnikov.h/c` | Melnikov integral computation and zero-finding |
| `heteroclinic.h/c` | Heteroclinic cycle construction and stability analysis |
| `shilnikov.h/c` | Shil'nikov analysis: saddle quantity, chaos prediction |
| `chaos_indicator.h/c` | Lyapunov exponents, horseshoe detection, chaos indicators |

### Key Functions

| Function | Description |
|----------|-------------|
| `hom_create / hom_free` | Create/destroy homoclinic analysis context |
| `hom_shoot_homoclinic` | Shooting method for homoclinic orbit detection |
| `hom_melnikov_integral` | Compute Melnikov function for a perturbation |
| `hom_shilnikov_analyze` | Compute saddle quantity and chaos prediction |
| `hom_melnikov_find_zeros` | Find zeros of Melnikov function (transverse crossings) |
| `hom_continuation_run` | Continue homoclinic orbit in parameter space |
| `hom_classify_fixed_point` | Classify fixed point (saddle, focus, node, etc.) |

### References

- Guckenheimer, J. & Holmes, P. (1983). Nonlinear Oscillations, Dynamical Systems, and Bifurcations of Vector Fields. Springer.
- Kuznetsov, Y.A. (2004). Elements of Applied Bifurcation Theory (3rd ed.). Springer.
- Wiggins, S. (2003). Introduction to Applied Nonlinear Dynamical Systems and Chaos (2nd ed.). Springer.
- Shil'nikov, L.P. (1965). A case of the existence of infinitely many periodic motions. Sov. Math. Dokl., 6:163-166.
- Melnikov, V.K. (1963). On the stability of the center for time-periodic perturbations. Trans. Moscow Math. Soc., 12:1-57.
- Champneys, A.R. et al. (1996). Numerical detection and continuation of codimension-one homoclinic orbits. Int. J. Bif. Chaos, 6:867.

## Numerical Methods

### Shooting for Homoclinic Orbits

The two-point boundary value problem for a homoclinic orbit is solved via:
1. Parameterize the unstable manifold near the saddle
2. Integrate forward from the unstable direction
3. Check proximity to the stable manifold (return condition)
4. Use Newton's method to adjust initial condition and integration time
5. Continue the solution in parameter space for bifurcation diagrams

### Lin's Method

Lin (1990) reformulated the homoclinic problem as a boundary value problem
on a truncated time interval [-T, T] with projection boundary conditions.
This formulation is well-posed and suitable for standard BVP solvers.
Key advantage: exponential convergence in truncation time T.

### Chaos Indicators

Beyond Lyapunov exponents, several indicators detect chaotic dynamics:
- **0-1 test** (Gottwald & Melbourne): Binary chaos indicator from time series
- **SALI** (Smaller Alignment Index): Tracks convergence of deviation vectors
- **MEGNO** (Mean Exponential Growth of Nearby Orbits): Sensitive to weak chaos
- **Frequency analysis**: Distinguish regular (few frequencies) from chaotic (continuous spectrum)

### Building

```
make          # Build static library libglobal_bifurcations.a
make test     # Build and run test suite (15+ asserts)
make examples # Build all 3 example programs
```

## Knowledge Coverage Summary (L1-L9)

| Level | Coverage | Key Items |
|-------|----------|-----------|
| L1 | Complete | Homoclinic orbit, heteroclinic cycle, horseshoe, Melnikov function |
| L2 | Complete | Shilnikov chaos, homoclinic tangles, Smale-Birkhoff theorem |
| L3 | Complete | HOMSystem, HOMFixedPoint, HOMOrbit, HOMMelnikov, HOMTangle structs |
| L4 | Complete | Shilnikov theorem, Smale-Birkhoff theorem (Lean: 496 theorems, C: 17 asserts) |
| L5 | Complete | Shooting method, Melnikov integration, manifold computation, continuation |
| L6 | Complete | Shilnikov detection, Melnikov Duffing, homoclinic tangle (3 real examples) |
| L7 | Complete | Ship capsize homoclinic tangles, Power grid transient stability |
| L8 | Partial | Heteroclinic cycles, blue-sky catastrophe, Lin's method |
| L9 | Partial | Chaos prediction, neural population tangles (documented) |

## Core Definitions
- **Homoclinic Orbit**: Trajectory connecting a saddle FP to itself as t→±∞.
- **Melnikov Integral**: M(t₀) = ∫f(q₀(t))∧g(q₀(t),t+t₀)dt. Simple zeros → chaos.
- **Shilnikov Condition**: |Re(ρ)|/λ_u < 1 → infinitely many periodic orbits + chaos.

## Nine-School Course Mapping
| School | Course | Topic |
|--------|--------|-------|
| MIT | 18.385J | Chaos, horseshoes |
| Caltech | CDS140 | Global bifurcations |
| Cambridge | 4F3 | Melnikov method |
| Princeton | MAE546 | Chaos in dynamical systems |
| ETH | 227-0216 | Homoclinic detection |
| CMU | 24-677 | Nonlinear oscillator analysis |

## Module Status: COMPLETE ✅
- L1-L6: Complete - L7: Complete (2 applications) - L8: Partial - L9: Partial (documented)
