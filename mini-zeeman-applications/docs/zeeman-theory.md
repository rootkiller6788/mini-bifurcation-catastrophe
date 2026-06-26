# Zeeman Catastrophe Theory — Core Framework

## Cusp Catastrophe

The cusp catastrophe is the most important of Thom's seven elementary catastrophes.
Its universal unfolding is:

```
V(x; a, b) = x^4/4 + a·x^2/2 + b·x
```

### Geometry

- **Equilibrium manifold**: M = {(x, a, b) | x^3 + a·x + b = 0}
- **Bifurcation set**: B = {(a, b) | 4a^3 + 27b^2 = 0}
- **Cusp point**: (a, b) = (0, 0)
- **Maxwell set**: b = 0 for a < 0 (equal-depth minima)

### Three Regions

| Region | Condition | Equilibria |
|--------|-----------|-----------|
| Outside cusp | Δ > 0 or a > 0 | 1 stable |
| Inside cusp | a < 0, Δ < 0 | 2 stable + 1 unstable |
| Bifurcation set | Δ = 0 | 1 stable + 1 degenerate |

### Zeeman's Five Qualities

1. **Bimodality**: Two distinct stable states coexist
2. **Sudden jumps**: Discontinuous state transitions
3. **Hysteresis**: Path-dependent behavior (forward ≠ backward)
4. **Divergence**: Small parameter changes → large state differences
5. **Inaccessibility**: Intermediate states are unstable

## Zeeman Catastrophe Machine

A rotating disk with elastic bands demonstrates all five qualities physically.

- Disk radius R, center at origin
- Fixed peg F, controllable peg C
- Elastic band from F and C to disk perimeter point P(θ)
- V(θ) = ½k₁|P-F|² + ½k₂|P-C|²
- As C moves, θ undergoes catastrophic jumps

## References

- Zeeman, E.C. (1972). "A Catastrophe Machine". Towards a Theoretical Biology, 4:276-282.
- Zeeman, E.C. (1976). "Catastrophe Theory". Scientific American, 234(4):65-83.
- Zeeman, E.C. (1977). Catastrophe Theory: Selected Papers 1972-1977. Addison-Wesley.
- Thom, R. (1975). Structural Stability and Morphogenesis. Benjamin.
- Poston, T. & Stewart, I. (1978). Catastrophe Theory and Its Applications. Pitman.
- Gilmore, R. (1981). Catastrophe Theory for Scientists and Engineers. Wiley.