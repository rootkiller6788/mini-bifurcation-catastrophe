# Thom Elementary Catastrophe Theory

## Rene Thom (1972) "Stabilite Structurelle et Morphogenese"

Thom proved that for gradient systems with <= 4 control parameters,
only seven structurally stable singularities exist:

### Cuspoids (1 state variable)
| Name | Potential V(x) | Codim |
|------|---------------|-------|
| Fold (A2) | x^3/3 + a*x | 1 |
| Cusp (A3) | x^4/4 + a*x^2/2 + b*x | 2 |
| Swallowtail (A4) | x^5/5 + a*x^3/3 + b*x^2/2 + c*x | 3 |
| Butterfly (A5) | x^6/6 + a*x^4/4 + b*x^3/3 + c*x^2/2 + d*x | 4 |

### Umbilics (2 state variables)
| Name | Potential V(x,y) | Codim |
|------|------------------|-------|
| Hyperbolic (D+4) | x^3/3 + y^3/3 + a*x*y + b*x + c*y | 3 |
| Elliptic (D-4) | x^3/3 - x*y^2 + a*(x^2+y^2) + b*x + c*y | 3 |
| Parabolic (D5) | x^2*y + y^4/4 + a*x^2 + b*y^2 + c*x + d*y | 4 |

## Key Concepts

### Equilibrium Manifold M
Set of points where gradient = 0: dV/dx = 0

### Bifurcation Set B
Projection of degenerate critical points (Hessian = 0) onto control space.
Divides control space into regions with different numbers of equilibria.

### Maxwell Set
Locus where two minima have equal depth. Relevant for hysteresis.

### Delay Convention
System stays in local minimum until it disappears, then jumps catastrophically.

### Maxwell Convention
System always occupies the global minimum.

## Zeeman Applications
- Heartbeat: cusp catastrophe model
- Brain: neural impulse propagation
- Ship stability: butterfly catastrophe
- Buckling: fold/cusp model
- Stock market crashes: cusp model
- Prison riots: social behavior cusp

## References
1. Thom (1972) Stabilite Structurelle et Morphogenese
2. Zeeman (1977) Catastrophe Theory: Selected Papers
3. Poston & Stewart (1978) Catastrophe Theory and Applications
4. Gilmore (1981) Catastrophe Theory for Scientists and Engineers
