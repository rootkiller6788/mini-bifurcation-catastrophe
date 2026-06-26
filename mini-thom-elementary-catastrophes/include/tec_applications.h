#ifndef TEC_APPLICATIONS_H
#define TEC_APPLICATIONS_H
#include "tec_core.h"

/* ==============================================================
 * tec_applications.h - Zeeman's Applications of Catastrophe Theory
 *
 * Zeeman applied catastrophe theory to diverse phenomena:
 *   1. Heartbeat: cusp catastrophe model of cardiac cycle
 *   2. Brain: cusp model of neural impulse propagation
 *   3. Prison riots: cusp model of social behavior
 *   4. Ship stability: butterfly catastrophe
 *   5. Buckling: fold/cusp model of Euler strut
 *   6. Stock market crashes: cusp with asymmetry
 *   7. Ecological regime shifts: fold bifurcation
 *   8. Phase transitions: Landau theory as cusp
 *
 * References:
 *   Zeeman (1977) Catastrophe Theory, Addison-Wesley
 *   Gilmore (1981) Catastrophe Theory, Wiley
 * ============================================================== */

/* --- Zeeman Heartbeat Model --- */
typedef struct { double diastole; double systole; double tension; double fiber_length; } TEC_HeartbeatState;
TEC_HeartbeatState tec_heartbeat_init(void);
void tec_heartbeat_step(TEC_HeartbeatState* hb, double dt);
double tec_heartbeat_get_ecg(const TEC_HeartbeatState* hb);

/* --- Euler Buckling Model (Fold Catastrophe) --- */
double tec_buckling_critical_load(double E, double I, double L);
double tec_buckling_deflection(double P, double P_crit, double imperfection);

/* --- Ship Stability (Butterfly Catastrophe) --- */
double tec_ship_restoring_moment(double heel_angle, double a, double b, double c, double d);
double tec_ship_capsize_angle(const double* butterfly_params);

/* --- Ecological Regime Shift (Fold-Cusp) --- */
double tec_regime_shift_threshold(double carrying_capacity, double harvesting_rate, double growth_rate);

/* --- Phase Transition (Landau) --- */
double tec_landau_free_energy(double order_param, double temperature, double Tc, double a4, double a6);
double tec_landau_equilibrium(double temperature, double Tc, double a4, double a6);

/* --- Stock Market Crash Model --- */
double tec_market_crash_probability(double sentiment, double leverage, double liquidity);
bool tec_market_is_in_crash_zone(double sentiment, double volatility, double leverage);

#endif
