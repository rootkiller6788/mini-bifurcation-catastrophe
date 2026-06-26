# Numerical Methods for Bifurcation Analysis

## Fixed Point Finding (Newton's Method)
x_{k+1} = x_k - J^{-1}*f(x_k). Quadratic convergence near root.
Requires good initial guess. Jacobian via finite differences.

## Continuation Methods
Natural parameter continuation: vary r, track fixed points.
Pseudo-arclength: parameterize by arclength s to pass folds.
Test functions: det(J) for SN, trace for Hopf.

## Eigenvalue Computation
2x2 analytic formula: lambda = (tr +/- sqrt(tr^2-4*det))/2.
For larger systems: numerical Jacobian + QR algorithm.

## Bifurcation Detection
Monitor eigenvalues during continuation.
Bracket root of test function, refine with secant method.
Check transversality and nondegeneracy conditions.

## Normal Form Coefficients
Saddle-node: a2 = (1/2)*d2f/dx2 (quadratic coefficient).
Hopf: l1 from Kuznetsov formula (first Lyapunov coefficient).
Pitchfork: a3 = (1/6)*d3f/dx3 (cubic coefficient).

## Implementation Architecture

### Layer 1: Linear Algebra (gst_core.c)
Vector/matrix operations, Gaussian elimination, eigenvalue computation.
Hand-written, no external dependencies. All edge cases handled.

### Layer 2: ODE Integration and Root Finding (gst_core.c)
RK4 for time integration. Newton-Raphson for fixed points.
Numerical Jacobian via finite differences.

### Layer 3: Bifurcation Detection (saddle_node.c, hopf_bifurcation.c, pitchfork_bifurcation.c)
Type-specific detection algorithms. Each uses eigenvalue crossing
patterns and test functions to identify bifurcation points.

### Layer 4: Analysis Tools
Normal form computation, continuation, stability analysis,
bifurcation diagrams, codimension-two analysis.
These compose the lower layers into complete analysis workflows.

## Performance Considerations
- 2x2 systems: analytic eigenvalues, O(1) per Jacobian
- NxN systems: O(N^3) for matrix inversion, O(N^2) for Jacobian
- Newton typically converges in < 10 iterations for good initial guesses
- Bifurcation detection requires ~50-200 parameter values
- Bifurcation diagram: ~200 parameter values, 5-10 fixed points each
