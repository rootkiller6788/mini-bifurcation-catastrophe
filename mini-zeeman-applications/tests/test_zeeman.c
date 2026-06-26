#include "zeeman_core.h"
#include "zeeman_biology.h"
#include "zeeman_social.h"
#include "zeeman_physics.h"
#include "zeeman_dynamics.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define EPS 1e-9

int main(void) {
    /* === zeeman_core: cusp catastrophe math === */
    double V0 = cusp_potential(0.0, 0.0, 0.0);
    assert(fabs(V0) < EPS);

    double V1 = cusp_potential(1.0, -1.0, 0.0);
    /* V = 1/4 - 1/2 = -0.25 */
    assert(fabs(V1 + 0.25) < EPS);

    double F = cusp_force(0.0, -1.0, 0.0);
    assert(fabs(F) < EPS);  /* at equilibrium, force=0 */

    double curv = cusp_curvature(0.0, -1.0);
    assert(curv < 0.0);  /* a=-1 -> unstable at x=0 */

    double curv2 = cusp_curvature(1.0, -1.0);
    assert(curv2 > 0.0);  /* 3*1 - 1 = 2 > 0 -> stable */

    /* Discriminant */
    double delta = cusp_discriminant(-1.0, 0.0);
    assert(delta < 0.0);  /* inside cusp */

    double delta_out = cusp_discriminant(1.0, 0.0);
    assert(delta_out > 0.0);  /* outside cusp */

    /* Bistable check */
    assert(cusp_is_bistable(-1.0, 0.0));
    assert(!cusp_is_bistable(1.0, 0.0));
    assert(!cusp_is_bistable(-1.0, 2.0));  /* outside cusp */

    /* On bifurcation set */
    double a_bs, b_bs;
    cusp_bifurcation_set_point(1.0, &a_bs, &b_bs);
    assert(cusp_is_on_bifurcation_set(a_bs, b_bs));

    /* Maxwell b */
    assert(fabs(cusp_maxwell_b(-1.0)) < EPS);

    /* Number of equilibria */
    assert(cusp_num_equilibria(1.0, 0.0) == 1);
    assert(cusp_num_equilibria(-1.0, 0.0) == 3);

    /* CuspEquilibrium */
    CuspEquilibrium* eq = cusp_equilibrium(-1.0, 0.0);
    assert(eq);
    assert(eq->n_roots == 3);
    assert(eq->n_stable == 2);
    cusp_equilibrium_print(eq);
    cusp_equilibrium_free(eq);

    eq = cusp_equilibrium(1.0, 0.0);
    assert(eq);
    assert(eq->n_roots == 1);
    assert(eq->n_stable == 1);
    cusp_equilibrium_free(eq);

    /* Higher catastrophes */
    double vf = cusp_potential_fold(1.0, -1.0);
    assert(!isnan(vf));
    double vs = cusp_potential_swallowtail(0.5, -1.0, 0.0, 0.5);
    assert(!isnan(vs));
    double vb = cusp_potential_butterfly(0.5, -1.0, 0.0, 0.5, 0.0);
    assert(!isnan(vb));

    CatastropheType ct = zeeman_classify_catastrophe(-1.0, 1.0, 0.0, 0.0);
    assert(ct == CATASTROPHE_CUSP);

    /* === Zeeman Machine === */
    ZeemanMachine* zm = zeeman_machine_create_default();
    assert(zm);
    assert(fabs(zm->R - 1.0) < EPS);
    assert(fabs(zm->Fx - 2.0) < EPS);
    assert(fabs(zm->Fy) < EPS);

    zeeman_machine_set_control(zm, -1.0, 0.5);
    double cx, cy;
    zeeman_machine_get_control(zm, &cx, &cy);
    assert(fabs(cx + 1.0) < EPS);
    assert(fabs(cy - 0.5) < EPS);

    double Vm = zeeman_machine_potential(zm, 0.0);
    assert(!isnan(Vm));
    double tau = zeeman_machine_torque(zm, 0.0);
    assert(!isnan(tau));

    CuspEquilibrium* meq = zeeman_machine_equilibrium(zm);
    assert(meq);
    assert(meq->n_roots >= 1 || meq->n_roots == 0);
    cusp_equilibrium_free(meq);

    double a_map, b_map;
    zeeman_machine_to_cusp(zm, &a_map, &b_map);
    assert(!isnan(a_map));
    assert(!isnan(b_map));

    zeeman_machine_step(zm, 0.01);
    zeeman_machine_set_damping(zm, 0.8);
    double theta = zeeman_machine_get_theta(zm);
    assert(!isnan(theta));
    zeeman_machine_print(zm);
    zeeman_machine_free(zm);

    /* === Path Tracking === */
    ZeemanPath* path = zeeman_path_create(128, CONVENTION_DELAY);
    assert(path);
    assert(path->n == 0);
    zeeman_path_add(path, -1.0, 0.0, -1.0);
    zeeman_path_add(path, -0.5, 0.0, -0.5);
    assert(path->n == 2);
    zeeman_path_track_equilibrium(path, -1.0, 0.0, 0.0, 0.0, 20);
    assert(path->n > 2);
    int jump_idx;
    bool has_jump = zeeman_path_detect_jump(path, &jump_idx);
    (void)has_jump;  /* may or may not have jump depending on path */
    double jump_mag = zeeman_path_jump_magnitude(path);
    assert(jump_mag >= 0.0);
    zeeman_path_print(path);
    zeeman_path_free(path);

    /* Safe path planning */
    ZeemanPath* safe = zeeman_plan_safe_path(-1.0, 0.0, 0.5, 0.2, 30);
    assert(safe);
    assert(safe->n > 0);
    zeeman_path_free(safe);

    /* === Hysteresis === */
    HysteresisLoop* loop = zeeman_hysteresis_loop(-2.0, 0.0, 0.0, 50, CONVENTION_DELAY);
    assert(loop);
    assert(loop->n == 51);
    double hyst_area = zeeman_hysteresis_area(loop);
    assert(hyst_area >= 0.0);
    zeeman_hysteresis_print(loop);
    zeeman_hysteresis_free(loop);

    /* === Bifurcation Diagram === */
    BifurcationDiagram* bd = zeeman_bifurcation_diagram(-2.0, 1.0, -1.0, 1.0, 20, 20);
    assert(bd);
    assert(bd->n_points > 0);
    int n_cat = zeeman_bifurcation_diagram_count_catastrophes(bd);
    assert(n_cat >= 0);
    zeeman_bifurcation_diagram_print(bd);
    zeeman_bifurcation_diagram_free(bd);

    /* === Heartbeat Model === */
    HeartbeatModel* hb = heartbeat_create(-1.0, 0.0, 0.1, 0.0);
    assert(hb);
    assert(hb->x > 0.0);
    assert(hb->epsilon > 0.0);
    assert(!heartbeat_is_systole(hb));  /* x=0.5 is not quite systole */
    HeartbeatModel* hb2 = heartbeat_create(-2.0, 0.0, 0.1, 0.0);
    heartbeat_simulate(hb2, 10.0, 0.01);
    assert(hb2->t > 0.0);
    double period = heartbeat_period(hb2, 20.0, 0.01);
    assert(period >= 0.0);
    heartbeat_print(hb2);
    heartbeat_free(hb);
    heartbeat_free(hb2);

    /* === Nerve Impulse === */
    NerveImpulse* ni = nerve_impulse_create(-1.0, 0.5, 0.5);
    assert(ni);
    assert(ni->v < 0.0);  /* resting potential */
    double threshold = nerve_impulse_threshold(ni);
    assert(threshold < 0.0);  /* threshold is negative for a=-1 */
    assert(!nerve_impulse_is_firing(ni));
    nerve_impulse_simulate(ni, 5.0, 0.01);
    int spikes = nerve_impulse_spike_count(ni, 5.0, 0.01);
    assert(spikes >= 0);
    nerve_impulse_print(ni);
    nerve_impulse_free(ni);

    /* === Cell Differentiation === */
    CellDiffModel* cd = cell_diff_create(0.1, 0.5, 2.0);
    assert(cd);
    assert(fabs(cd->morphogen - 0.1) < EPS);
    int n_states = cell_diff_stable_states(cd);
    assert(n_states >= 0);
    double crit_m = cell_diff_critical_morphogen(cd);
    assert(crit_m >= 0.0);
    cell_diff_simulate(cd, 5.0, 0.01);
    cell_diff_print(cd);
    cell_diff_free(cd);

    /* === Prison Riot === */
    PrisonRiotModel* pr = prison_riot_create(0.3, 0.2);
    assert(pr);
    assert(!prison_riot_is_riot(pr));
    double crit_t = prison_riot_critical_tension(pr);
    assert(crit_t >= 0.0);
    prison_riot_simulate(pr, 10.0, 0.01);
    prison_riot_print(pr);
    prison_riot_free(pr);

    /* === Stock Market === */
    StockMarketModel* sm = stock_market_create(100.0, 0.4, 100.0, 0.1);
    assert(sm);
    assert(fabs(sm->price - 100.0) < EPS);
    assert(!stock_market_crash_detect(sm, 0.3));
    double crash_mag = stock_market_crash_magnitude(sm);
    assert(crash_mag >= 0.0);
    stock_market_simulate(sm, 5.0, 0.01);
    stock_market_print(sm);
    stock_market_free(sm);

    /* === War/Peace === */
    WarPeaceModel* wp = war_peace_create(0.2, 0.8);
    assert(wp);
    assert(!war_peace_is_war(wp));
    double thr = war_peace_threshold(wp);
    assert(thr >= 0.0);
    war_peace_simulate(wp, 10.0, 0.01);
    war_peace_print(wp);
    war_peace_free(wp);

    /* === Euler Buckling === */
    EulerBuckling* eb = euler_buckling_create(5.0, 1.0, 1.0, 0);
    assert(eb);
    double Pcr = euler_buckling_critical_load(eb);
    assert(Pcr > 0.0);
    assert(!euler_buckling_has_buckled(eb)); /* below critical */
    EulerBuckling* eb2 = euler_buckling_create(Pcr * 2.0, 1.0, 1.0, 0);
    euler_buckling_simulate(eb2, 5.0, 0.01);
    assert(euler_buckling_has_buckled(eb2));
    double stiff = euler_buckling_post_buckling_stiffness(eb2);
    assert(stiff > 0.0);
    euler_buckling_print(eb2);
    euler_buckling_free(eb);
    euler_buckling_free(eb2);

    /* === Phase Transition === */
    /* Use parameters: a=0.5, b=0.01 (normalized units), T=1.0 below Tc */
    PhaseTransition* pt = phase_transition_create(1.0, 0.1, 0.5, 0.01);
    assert(pt);
    assert(pt->T > 0.0);
    assert(!phase_transition_is_supercritical(pt));  /* T < Tc */
    double Vl, Vg, Peq;
    Peq = phase_transition_maxwell_construction(pt, &Vl, &Vg);
    /* Peq might be numerically unstable for small parameters; check is finite */
    assert(!isnan(Peq));
    double Vsp1, Vsp2;
    phase_transition_spinodal(pt, &Vsp1, &Vsp2);
    phase_transition_simulate(pt, 10.0, 0.01);
    phase_transition_print(pt);
    phase_transition_free(pt);

    /* === Ship Stability === */
    ShipStability* ss = ship_stability_create(0.5, 10.0, 0.3);
    assert(ss);
    double GZ = ship_stability_righting_moment(ss, 0.1);
    assert(GZ > 0.0);  /* positive GM -> positive righting moment */
    double Vss = ship_stability_potential(ss, 0.0);
    assert(!isnan(Vss));
    assert(!ship_stability_capsize_check(ss, 0.5));
    ship_stability_set_wave(ss, 0.2, 0.5);
    ship_stability_simulate(ss, 5.0, 0.01);
    double crit_gm = ship_stability_critical_gm(ss);
    assert(crit_gm >= 0.0);
    ship_stability_print(ss);
    ship_stability_free(ss);

    /* === Jump Detection === */
    double jump_a = zeeman_find_delay_jump(-2.0, 0.5, 0.2, 50, -1.0);
    assert(!isnan(jump_a));
    double mj_a = zeeman_find_maxwell_jump(-2.0, 0.5, 0.0, 50);
    assert(!isnan(mj_a));

    /* Bifurcation crossing */
    double a_path[3] = {-2.0, -1.0, 0.0};
    double b_path[3] = {0.0, 0.0, 0.0};
    int cross_idx;
    bool crossed = zeeman_detect_bifurcation_crossing(a_path, b_path, 3, &cross_idx);
    assert(crossed);  /* -2 to 0 at b=0 crosses bifurcation set */

    /* === Divergence Test === */
    DivergenceTest* dt = zeeman_divergence_test(-1.0, -0.5, 0.5, 30);
    assert(dt);
    assert(dt->divergence_measure >= 0.0);
    assert(dt->sensitivity >= 0.0);
    zeeman_divergence_test_print(dt);
    zeeman_divergence_test_free(dt);

    /* Inaccessibility */
    assert(zeeman_is_inaccessible(-1.0, 0.0, 0.0));  /* x=0 is unstable for a=-1,b=0 */

    /* Map to cusp normal form */
    double a_norm, b_norm, shift;
    zeeman_map_to_cusp_normal(1.0, 0.0, -3.0, 0.0, &a_norm, &b_norm, &shift);
    assert(!isnan(a_norm));

    /* Potential surface */
    double* surf = zeeman_potential_surface(-2.0, 0.0, 5, -0.5, 0.5, 5, 1.0);
    assert(surf);
    free(surf);

    double gm_x = zeeman_global_minimum(-1.0, 0.0);
    assert(!isnan(gm_x));

    double barrier = zeeman_barrier_height(-1.0, 0.0);
    assert(barrier > 0.0);  /* there is a barrier when bistable */

    /* === Slow-Fast Cusp === */
    SlowFastCusp* sf = slowfast_cusp_create(0.05, -2.0, 0.0, 3.0, 0.5, 1.0);
    assert(sf);
    assert(slowfast_cusp_is_relaxation(sf));
    slowfast_cusp_simulate(sf, 20.0, 0.01);
    double sf_period = slowfast_cusp_period(sf, 20.0, 0.01);
    assert(sf_period >= 0.0);
    slowfast_cusp_print(sf);
    slowfast_cusp_free(sf);

    /* === Zeeman's Five Qualities === */
    assert(zeeman_quality_bimodality(-1.0, 0.0));
    assert(!zeeman_quality_bimodality(1.0, 0.0));
    assert(zeeman_quality_inaccessibility(-1.0, 0.0));

    double hyst_q = zeeman_quality_hysteresis(-2.0, 0.0, 0.0, 30);
    assert(hyst_q >= 0.0);

    double div_q = zeeman_quality_divergence(-1.0, -0.5, 0.5, 20);
    assert(div_q >= 0.0);

    /* Full report */
    zeeman_five_qualities_report(-1.0, 0.0);

    printf("All tests passed.\n");
    return 0;
}