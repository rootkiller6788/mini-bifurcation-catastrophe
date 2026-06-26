# Homoclinic Orbits and Global Bifurcations

## Homoclinic Orbits
A homoclinic orbit is a trajectory that connects a saddle fixed point to itself as t->inf and t->-inf. It lies in the intersection of the stable and unstable manifolds.

## Shilnikov Bifurcation (1965)
When a homoclinic orbit exists at a saddle-focus (one real eigenvalue, two complex), the Shilnikov condition determines whether chaos emerges:
- |Re(lambda_s)| < Re(lambda_u): infinitely many periodic orbits, chaos
- |Re(lambda_s)| > Re(lambda_u): single stable periodic orbit

## Melnikov Method (1963)
Analytical perturbation method for detecting homoclinic bifurcations. The Melnikov integral M(t0) measures the splitting distance between stable and unstable manifolds. Simple zeros of M(t0) indicate transverse homoclinic intersections.

## Smale Horseshoe
Transverse homoclinic intersections imply the existence of a Smale horseshoe map, which has:
- A Cantor set invariant set
- Dense periodic orbits
- Sensitive dependence on initial conditions
- Positive topological entropy

## Global Bifurcations
Unlike local bifurcations (detected via eigenvalues), global bifurcations involve changes in the global phase portrait structure. They include homoclinic/heteroclinic bifurcations, blue sky catastrophes, and homoclinic explosions.
