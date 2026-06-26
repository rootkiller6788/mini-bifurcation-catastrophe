# Delay Differential Equations: Theory

## Retarded Functional Differential Equations
dx/dt = f(x(t), x(t-tau), mu)
The state derivative depends on both current and delayed states.
Requires history function on [-tau, 0] for well-posedness.

## Method of Steps
Integrate ODE on intervals [0,tau], [tau,2tau], ...
Solution is pieced together from ODE segments.
Numerical: store history in circular buffer.

## Characteristic Equation
Linearized DDE: dx/dt = A*x(t) + B*x(t-tau)
Characteristic: det(lambda*I - A - B*exp(-lambda*tau)) = 0
Quasipolynomial with infinitely many roots.

## Lambert W Function
Solves w*exp(w) = z
Principal branch W_0: real-valued for z >= -1/e
Used for analytical stability of scalar DDEs.

## Stability Switches
As tau varies, roots cross imaginary axis.
Stability lobes: alternating stable/unstable regions.
Number of unstable roots changes only at crossings.

## Hopf Bifurcation in DDEs
Pair of complex conjugate roots crosses imaginary axis.
Normal form yields Lyapunov coefficient.
Supercritical: stable limit cycle born at bifurcation.

## Delayed Feedback Control
u(t) = K*(x(t-tau_c) - x(t))
Noninvasive: vanishes on tau_c-periodic orbit.
Stability determined by Floquet multipliers of controlled system.

## Amplitude Death
Diffusive coupling with delay quenches oscillations.
Critical coupling strength depends on delay and topology.
Hysteresis: oscillations restart at lower coupling than death.
