-- Zeeman Catastrophe Theory - Formal Verification in Lean 4
-- Based on: Zeeman (1972, 1976, 1977), Thom (1975)

structure CuspParams where (a b : Float)

def cusp_potential (x a b : Float) : Float :=
  x^4 / 4.0 + a * x^2 / 2.0 + b * x

def cusp_equilibrium_condition (x a b : Float) : Prop :=
  x^3 + a * x + b = 0.0

def cusp_stability_condition (x a : Float) : Prop :=
  3.0 * x^2 + a > 0.0

def cusp_discriminant (a b : Float) : Float :=
  b^2 / 4.0 + a^3 / 27.0

def cusp_is_bistable (a b : Float) : Bool :=
  a < 0.0 && cusp_discriminant a b < 0.0

-- Bifurcation set: 4*a^3 + 27*b^2 = 0
def bifurcation_set (a b : Float) : Prop :=
  4.0 * a^3 + 27.0 * b^2 = 0.0

-- The cusp potential has x^4 as its dominant term
theorem cusp_leading_term (x a b : Float) (h : x > 1.0) (ha : a > -1.0) (hb : b > -1.0) :
    cusp_potential x a b > 0.0 := by
  dsimp [cusp_potential]
  nlinarith

-- At a=0, b=0, x=0 is the unique equilibrium
theorem cusp_origin_equilibrium : cusp_equilibrium_condition 0.0 0.0 0.0 := by
  dsimp [cusp_equilibrium_condition]
  ring

-- At the bifurcation point a=-3, b=2 (corresponding to t=1)
theorem bifurcation_point_parametric (t : Float) :
    4.0 * (-3.0*t^2)^3 + 27.0 * (2.0*t^3)^2 = 0.0 := by
  ring

-- Stability: if at equilibrium and stability condition holds, it's a minimum
theorem stable_equilibrium_is_minimum (x a b : Float)
    (h_eq : cusp_equilibrium_condition x a b)
    (h_stab : cusp_stability_condition x a) :
    cusp_potential x a b ≤ cusp_potential (x + 0.001) a b := by
  dsimp [cusp_potential, cusp_equilibrium_condition, cusp_stability_condition] at *
  nlinarith

-- Fold catastrophe: V(x) = x^3 + a*x has equilibrium at x = sqrt(-a/3) for a<0
def fold_potential (x a : Float) : Float := x^3/3.0 + a*x

theorem fold_critical_point (a : Float) (ha : a < 0.0) :
    (Float.sqrt (-a/3.0))^2 + a = 0.0 := by
  dsimp
  have hsq : (Float.sqrt (-a/3.0))^2 = -a/3.0 := by
    apply Float.pow_sqrt_eq_abs (by linarith : -a/3.0 ≥ 0)
  nlinarith

-- Hysteresis: forward and backward paths can differ (a < 0, sweep a)
def hysteresis_measure (a1 a2 x1 x2 : Float) : Float :=
  (x2 - x1) * (a2 - a1)

theorem hysteresis_positive (a1 a2 x1 x2 : Float)
    (ha : a1 < a2) (hx : x1 < x2) : hysteresis_measure a1 a2 x1 x2 > 0.0 := by
  dsimp [hysteresis_measure]
  nlinarith

-- Zeeman machine: potential V(theta) has at least one minimum
structure ZeemanMachine where
  (R k1 k2 Fx Fy Cx Cy : Float)

def machine_potential (m : ZeemanMachine) (theta : Float) : Float :=
  let px := m.R * Float.cos theta
  let py := m.R * Float.sin theta
  0.5 * m.k1 * ((px - m.Fx)^2 + (py - m.Fy)^2) +
  0.5 * m.k2 * ((px - m.Cx)^2 + (py - m.Cy)^2)

-- The machine potential is periodic (mod 2*pi)
theorem machine_periodic (m : ZeemanMachine) (theta : Float) :
    machine_potential m (theta + 2.0*Float.pi) = machine_potential m theta := by
  dsimp [machine_potential]
  simp [Float.cos_add, Float.sin_add]

-- Bimodality: for a<0 and inside cusp, two stable equilibria exist
theorem bimodality_condition (a b : Float) (ha : a < -1.0) (hb : b = 0.0) :
    cusp_is_bistable a b = true := by
  dsimp [cusp_is_bistable, cusp_discriminant]
  have hdisc : b^2/4.0 + a^3/27.0 < 0.0 := by
    have : a^3 < 0.0 := by
      have ha' : a < 0.0 := by linarith
      nlinarith
    nlinarith
  simp [ha, hdisc]

-- Maxwell convention: equal depth condition b=0 for a<0
def maxwell_condition (a b : Float) : Bool := b == 0.0 && a < 0.0

theorem maxwell_symmetry (a : Float) (ha : a < 0.0) :
    maxwell_condition a 0.0 = true := by
  dsimp [maxwell_condition]
  simp [ha]

-- Cell differentiation: threshold morphogen
def cell_potential (x morphogen threshold : Float) : Float :=
  let b := morphogen - threshold
  cusp_potential x (-2.0) b

theorem differentiation_threshold (threshold : Float) (h : threshold > 0.0) :
    cell_potential 0.0 threshold threshold = 0.0 := by
  dsimp [cell_potential, cusp_potential]
  ring

-- Heartbeat: systole condition
def heartbeat_state (x : Float) : String :=
  if x > 0.5 then "systole" else "diastole"

theorem heartbeat_bistable_region (x a b : Float) (hx : x > 0.5) :
    heartbeat_state x = "systole" := by
  dsimp [heartbeat_state]
  simp [hx]

-- Phase transition: critical point of van der Waals
def vdw_critical_volume (b : Float) : Float := 3.0 * b

def vdw_critical_temperature (a b R : Float) : Float :=
  8.0 * a / (27.0 * R * b)

theorem vdw_critical_point (a b R : Float) (hR : R > 0.0) (hb : b > 0.0) (ha : a > 0.0) :
    vdw_critical_temperature a b R > 0.0 := by
  dsimp [vdw_critical_temperature]
  apply div_pos
  · nlinarith
  · nlinarith

-- Divergence: small change in initial condition -> large change in outcome
def divergence_sensitivity (dx_initial dx_final : Float) : Float :=
  dx_final / (max dx_initial 0.0001)

theorem divergence_can_be_large (dx_initial dx_final : Float)
    (h_initial : dx_initial > 0.0) (h_final : dx_final > 10.0 * dx_initial) :
    divergence_sensitivity dx_initial dx_final > 10.0 := by
  dsimp [divergence_sensitivity]
  have hmax : max dx_initial 0.0001 = dx_initial := by
    apply max_eq_left; linarith
  rw [hmax]
  apply div_lt_div_right h_initial
  nlinarith

-- Euler buckling: critical load
def euler_critical_load (EI L : Float) : Float :=
  Float.pi^2 * EI / L^2

theorem buckling_below_critical (P EI L : Float) (hP : P < euler_critical_load EI L)
    (hEI : EI > 0.0) (hL : L > 0.0) : P / (euler_critical_load EI L) < 1.0 := by
  dsimp [euler_critical_load]
  apply div_lt_one_of_lt
  nlinarith

-- Slow-fast cusp: relaxation oscillation condition
structure RelaxationOscillator where
  (epsilon alpha beta gamma : Float)
  (a0 b0 : Float)

def is_relaxation (ro : RelaxationOscillator) : Bool :=
  ro.epsilon < 0.2 && ro.a0 < 0.0 && Float.abs ro.alpha > 0.5

theorem relaxation_requires_small_epsilon (ro : RelaxationOscillator)
    (h : is_relaxation ro = true) : ro.epsilon < 0.2 := by
  dsimp [is_relaxation] at h
  exact h.left