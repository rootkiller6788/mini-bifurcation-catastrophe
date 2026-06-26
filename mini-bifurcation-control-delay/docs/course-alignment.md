# Course Alignment -- Bifurcation Control & Delay

## MIT 18.385J: Nonlinear Dynamics and Chaos
- **Delay differential equations**: dx/dt = f(x(t), x(t-tau))
- **Hopf bifurcation in DDE**: tau as bifurcation parameter

## Caltech CDS140: Nonlinear Dynamics
- **Washout filter**: Stabilize bifurcation without changing equilibrium
- **Delayed feedback control (Pyragas)**: u(t) = K[y(t) - y(t-tau)]

## Key Formulas
```
Washout filter: dw/dt = x - omega*w, u = -k(x - omega*w)
DFC: u(t) = K[y(t-tau) - y(t)]
Characteristic equation (DDE): det(lambda*I - A - B*exp(-lambda*tau)) = 0
```

## Exercises
1. Compute critical delay for Hopf bifurcation in DDE
2. Design washout filter to stabilize saddle-node bifurcation
3. Apply Pyragas DFC to stabilize UPO in Lorenz system
