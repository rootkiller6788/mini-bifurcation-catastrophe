# Course Alignment -- Global Bifurcations & Homoclinic Orbits

## MIT 18.385J: Nonlinear Dynamics and Chaos
- **Homoclinic orbit**: Trajectory connecting saddle to itself as t->+-inf
- **Shilnikov chaos**: Homoclinic orbit to saddle-focus -> Smale horseshoe

## Caltech CDS140: Nonlinear Dynamics
- **Melnikov method**: Perturbative detection of homoclinic tangles
- **Horseshoe map**: Symbolic dynamics proof of chaos

## Key Formulas
```
Melnikov function: M(t0) = integral_{-inf}^{inf} f0 ^ g(q_u0(t)) dt
  M(t0) has simple zero -> homoclinic tangle -> chaos

Shilnikov condition (3D): |Re(rho)| < Re(lambda_u)
  (saddle-focus eigenvalues: -rho+-i*omega, lambda_u>0)

Saddle quantity: sigma = -Re(rho) + lambda_u
  sigma<0 -> single orbit, sigma>0 -> infinite orbits
```

## Exercises
1. Compute Melnikov function for forced Duffing oscillator
2. Verify Shilnikov condition for Chua's circuit
3. Numerically trace homoclinic orbit continuation
