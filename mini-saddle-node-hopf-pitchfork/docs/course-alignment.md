# Course Alignment -- Saddle-Node, Hopf & Pitchfork Bifurcations

## MIT 18.385J: Nonlinear Dynamics and Chaos
- **Saddle-node (fold)**: Normal form dx/dt = r + x^2. Bifurcation at r=0.
- **Hopf bifurcation**: Normal form dz/dt = (alpha+i*omega)z + (a+i*b)z|z|^2. Supercritical if a<0.

## Caltech CDS140: Nonlinear Dynamics
- **Codimension**: Minimal number of parameters to unfold bifurcation
- **Bifurcation diagram**: Fixed points/limit cycles vs. parameter

## Key Formulas
```
Saddle-node: dx/dt = r + x^2
  r<0: stable + unstable FP . r=0: semi-stable . r>0: no FP

Supercritical Hopf: dz/dt = (alpha+i*omega)z - |z|^2*z  (a<0)
  alpha<0: stable FP . alpha=0: weak focus . alpha>0: unstable FP + stable limit cycle

Supercritical pitchfork: dx/dt = r*x - x^3
  r<0: one stable FP (x=0) . r>0: unstable FP (x=0) + two stable FPs (+-sqrt(r))
```

## Exercises
1. Detect saddle-node bifurcation in 1D system
2. Compute first Lyapunov coefficient for Hopf bifurcation
3. Generate complete bifurcation diagram
