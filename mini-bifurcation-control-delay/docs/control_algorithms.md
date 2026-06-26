# Bifurcation Control Algorithms for DDEs

## DDE Integration (RK4)
1. Maintain circular buffer of past states
2. For each step: retrieve x(t) from buffer, x(t-tau) by indexing
3. Evaluate f(x(t), x(t-tau), t) for RK4 slopes
4. Advance and store new state in buffer

## Hopf Detection Algorithm
1. Vary bifurcation parameter mu
2. For each mu: reset DDE, simulate transient
3. Linearize at equilibrium: compute Jacobians J0, J_tau
4. Solve characteristic equation numerically
5. Detect sign change in max(Re(lambda))

## Pyragas DFC Design
1. Estimate UPO period tau_c (from FFT or Poincare)
2. Set control delay tau_control = tau_c
3. Sweep gain K: for each K, test stability
4. Optimal K from Nyquist criterion or Floquet analysis

## Washout Filter Design
1. Identify oscillation frequency omega from FFT
2. Set washout pole omega_w ≈ omega
3. Choose alpha for desired bandwidth
4. Cascade with DFC for robust stabilization

## Stability Chart Computation
1. Grid (tau, K) parameter space
2. For each point: compute rightmost eigenvalue
3. Boundary where max(Re(lambda)) = 0 is stability boundary
4. Count stability lobes

## Amplitude Death Detection
1. Initialize coupled DDEs with oscillatory initial conditions
2. Simulate with increasing coupling strength
3. Monitor oscillation amplitude over time
4. Death detected when amplitude < threshold
5. Record critical coupling and death time
