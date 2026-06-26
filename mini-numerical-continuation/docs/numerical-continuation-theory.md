# Numerical Continuation Theory

## The Problem

Given F(x, lambda)=0 with known solution (x0, lambda0), trace the solution
branch as lambda varies. At turning points (folds), lambda is not a valid
continuation parameter.

## Pseudo-Arclength Continuation (Keller 1977)

Augment F(x,lambda)=0 with arclength constraint:
  (x - x_k) * t_x + (lambda - lambda_k) * t_lambda = ds

where (t_x, t_lambda) is the unit tangent at point k.

## Predictor-Corrector

1. Predictor: extrapolate along tangent x_pred = x_k + ds * t
2. Corrector: Newton iteration on augmented system until convergence

## Bifurcation Detection

Test functions along the branch:
- Limit point: sign change in d(lambda)/ds
- Branch point: sign change in det(J)
- Hopf: pair of eigenvalues crosses imaginary axis

## Step Size Control

Adapt ds based on Newton convergence:
- Fast convergence → increase ds
- Slow/no convergence → decrease ds
- Optimal: 4-6 Newton iterations per step

## References
- Keller (1977) Numerical Solution of Bifurcation Problems
- Allgower & Georg (2003) Numerical Continuation Methods
- Seydel (2010) Practical Bifurcation and Stability Analysis
