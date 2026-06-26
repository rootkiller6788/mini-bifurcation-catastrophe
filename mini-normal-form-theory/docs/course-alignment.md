# Course Alignment -- Normal Form Theory

## MIT 18.385J: Nonlinear Dynamics and Chaos
- **Poincare-Birkhoff normal form**: Simplify nonlinear system near equilibrium
- **Resonance condition**: m*lambda = 0 prevents removing certain nonlinear terms

## Caltech CDS140: Nonlinear Dynamics
- **Center manifold reduction**: Reduce dimension near bifurcation
- **Unfolding**: Adding parameters to capture all perturbations

## Key Formulas
```
Linear part:  dx/dt = Ax + f2(x) + f3(x) + ...
Transform:    x = y + h(y), h homogeneous polynomial
Homological equation: L_A(h_k) = f_k - g_k

Normal form for saddle-node: dx/dt = r +- x^2
Normal form for Hopf: dz/dt = (alpha+i*omega)z + (a+i*b)z|z|^2
```

## Exercises
1. Compute normal form for saddle-node bifurcation
2. Derive normal form for Hopf bifurcation (compute a,b coefficients)
3. Implement center manifold approximation
