# Course Alignment — Numerical Continuation

## MIT 6.241J: Dynamic Systems and Control
- Equilibrium tracking as parameters vary
- Bifurcation analysis of control systems

## Stanford AA203: Optimal and Learning-based Control
- Continuation for solving optimality conditions
- Homotopy methods for non-convex optimization

## Berkeley EE222: Nonlinear Systems
- Bifurcation diagrams via continuation
- Turning points and branch switching

## Caltech CDS110: Introduction to Control Theory
- Root locus as a continuation problem
- Gain margin from bifurcation analysis

## Key Formulas
- F(x,lambda)=0
- J*dx + F_lambda*dlambda = 0
- (x-x0)*t_x + (lambda-lambda0)*t_lambda = ds
- Newton: x_{k+1} = x_k - J^{-1}*F(x_k)
- Feigenbaum delta = lim(d_n/d_{n+1}) = 4.669
