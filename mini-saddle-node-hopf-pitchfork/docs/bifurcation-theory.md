# Local Bifurcation Theory

## Saddle-Node (Fold) Bifurcation
Normal form: dx/dt = r +/- x^2
Eigenvalue condition: simple zero eigenvalue crossing real axis.
r<0: two fixed points (node + saddle). r=0: one semi-stable FP. r>0: no FPs.
Codimension: 1. Unfolding parameter: r.

## Hopf Bifurcation
Normal form: dz/dt = (r+i*omega)*z + (a+i*b)*|z|^2*z
Eigenvalue condition: complex conjugate pair crosses imaginary axis.
First Lyapunov coefficient l1: l1<0 supercritical (stable LC), l1>0 subcritical.
Codimension: 1.

## Pitchfork Bifurcation
Normal form: dx/dt = r*x +/- x^3
Requires Z2 symmetry: f(-x,r) = -f(x,r).
Supercritical (x^3 term negative): continuous transition.
Subcritical (x^3 term positive): discontinuous with hysteresis.
Codimension: 1.

## Center Manifold Theorem
Near a bifurcation, the dynamics is topologically equivalent to
the restriction to the center manifold. This reduces dimensionality.

## References
Kuznetsov (2004), Strogatz (2015), Guckenheimer & Holmes (1983)

## Codimension-Two Bifurcations

### Bogdanov-Takens
Double zero eigenvalue. Normal form requires 2 parameters.
Produces SN, Hopf, and homoclinic bifurcation curves.

### Cusp
Two SN curves meet at a cusp. Organizing center for hysteresis.
Normal form: x' = r + sigma*x - x^3.

### Bautin (Generalized Hopf)
First Lyapunov coefficient vanishes. Separates super/subcritical Hopf.

### Fold-Hopf (Zero-Hopf)
Zero eigenvalue + purely imaginary pair. Produces quasiperiodicity.

## Numerical Methods

### Continuation
Pseudo-arclength method parameterizes solution curve by arclength s,
allowing continuation through folds where natural parameter continuation fails.

### Test Functions
det(J) for saddle-node detection. trace(J) for Hopf detection.
Monitor for sign changes, then refine with secant method.

### Branch Switching
At bifurcation points, switch to new solution branch using
eigenvector directions as initial guesses.
