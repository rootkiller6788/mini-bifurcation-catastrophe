# Bifurcation & Catastrophe Theory

A collection of **from-scratch, zero-dependency C implementations** of bifurcation and catastrophe theory for nonlinear dynamical systems. Each sub-module maps to MIT, Caltech, Stanford and other top-tier university courses, bridging theory and practice by translating textbook equations into runnable C code. Based on Kuznetsov (2004), Thom (1975), Zeeman (1977).

## Sub-Modules

| Sub-Module | Topics | Key Courses |
|--------|--------|-------------|
| [mini-saddle-node-hopf-pitchfork](mini-saddle-node-hopf-pitchfork/) | Saddle-node, Hopf, pitchfork bifurcations; center manifold reduction; codimension-2 analysis; eigenvalue crossing detection | MIT 18.385J, Caltech CDS140 |
| [mini-normal-form-theory](mini-normal-form-theory/) | Normal form computation, center manifold reduction (approximation), resonance analysis, universal unfolding | MIT 18.385J, Cornell MAE 5790 |
| [mini-numerical-continuation](mini-numerical-continuation/) | Pseudo-arclength continuation (Keller), predictor-corrector, branch switching, Moore-Penrose continuation | MIT 18.085, Caltech CDS140 |
| [mini-global-bifurcations-homoclinic](mini-global-bifurcations-homoclinic/) | Homoclinic/heteroclinic orbits, Shilnikov theorems, Melnikov method, Lin's method, figure-eight orbits | MIT 18.385J, Princeton MAE 546 |
| [mini-bifurcation-in-pde-patterns](mini-bifurcation-in-pde-patterns/) | Turing instability, reaction-diffusion PDE bifurcation, linear stability analysis, pattern formation | MIT 18.385J, Stanford MATH 220 |
| [mini-bifurcation-control-delay](mini-bifurcation-control-delay/) | Delay differential equations (DDE), Pyragas delayed feedback control (DFC), Lambert W stability, Hopf bifurcation in delay systems | MIT 6.832, ETH 227-0216 |
| [mini-thom-elementary-catastrophes](mini-thom-elementary-catastrophes/) | Thom's seven elementary catastrophes: fold, cusp, swallowtail, butterfly, hyperbolic/elliptic/parabolic umbilics | MIT 18.385J, Harvard Math 118 |
| [mini-zeeman-applications](mini-zeeman-applications/) | Zeeman catastrophe machine, cusp geometry, biological/social/economic catastrophe models | MIT 18.385J, Cambridge Part III |

## Design Philosophy

- **Zero external dependencies** — pure C (C99/C11), only `libc` and `libm`
- **Self-contained modules** — each directory has its own `Makefile`, `include/`, `src/`, `examples/`, `demos/`, `tests/`
- **Theory-to-code mapping** — every sub-module includes `docs/` with course-alignment notes
- **Practical demos** — bifurcation detectors, continuation solvers, catastrophe surface visualizers, and more

## Building

Each sub-module is standalone. Navigate to a sub-module directory and run:

```bash
cd mini-bifurcation-control-delay
make all    # build everything
make test   # run tests
```

Requires **GCC** and **GNU Make**.

## Project Structure

```
mini-bifurcation-catastrophe/
├── mini-saddle-node-hopf-pitchfork/  # Saddle-node, Hopf, pitchfork bifurcations; center manifold
├── mini-normal-form-theory/          # Normal form computation and universal unfolding
├── mini-numerical-continuation/      # Numerical continuation and branch tracing
├── mini-global-bifurcations-homoclinic/ # Global bifurcations: homoclinic orbits, Melnikov method
├── mini-bifurcation-in-pde-patterns/ # Bifurcation analysis in PDE systems and pattern formation
├── mini-bifurcation-control-delay/   # Bifurcation control under delay differential equations
├── mini-thom-elementary-catastrophes/ # Thom's seven elementary catastrophes
└── mini-zeeman-applications/         # Zeeman catastrophe machine and real-world applications
```

## License

MIT
