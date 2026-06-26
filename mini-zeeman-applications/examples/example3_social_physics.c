#include "zeeman_social.h"
#include "zeeman_physics.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=======================================\n");
    printf("  Zeeman Social & Physics — Demo 3\n");
    printf("  Prison Riots, Stock Market, Buckling\n");
    printf("=======================================\n\n");

    /* 1. Prison riot */
    printf("--- 1. Prison Riot Model ---\n");
    PrisonRiotModel* pr_low = prison_riot_create(0.2, 0.5);
    prison_riot_simulate(pr_low, 20.0, 0.01);
    printf("Low tension (0.2): "); prison_riot_print(pr_low);
    prison_riot_free(pr_low);

    PrisonRiotModel* pr_high = prison_riot_create(0.8, 0.7);
    prison_riot_simulate(pr_high, 20.0, 0.01);
    printf("High tension (0.8): "); prison_riot_print(pr_high);
    prison_riot_free(pr_high);

    /* 2. Stock market */
    printf("\n--- 2. Stock Market Crash Model ---\n");
    StockMarketModel* sm_bull = stock_market_create(120.0, 0.7, 100.0, 0.15);
    printf("Bull market: "); stock_market_print(sm_bull);
    printf("Simulating 30 time units with high speculation...\n");
    for (int i = 0; i < 30; i++) {
        stock_market_step(sm_bull, 1.0);
        if (i % 10 == 0) {
            printf("  t=%.0f: price=%.2f spec=%.3f crash=%s\n",
                   sm_bull->t, sm_bull->price, sm_bull->speculation,
                   stock_market_crash_detect(sm_bull, 0.2) ? "DETECTED" : "none");
        }
    }
    printf("Final: "); stock_market_print(sm_bull);
    stock_market_free(sm_bull);

    /* 3. War/peace */
    printf("\n--- 3. War/Peace Model ---\n");
    printf("Escalating threat from 0.3 to 0.7:\n");
    WarPeaceModel* wp = war_peace_create(0.3, 0.5);
    for (int i = 0; i <= 10; i++) {
        wp->threat = 0.3 + 0.4 * i / 10.0;
        wp->a_factor = -(wp->threat - 0.5);
        war_peace_step(wp, 1.0);
        printf("  threat=%.2f state=%.3f war=%s\n",
               wp->threat, wp->state,
               war_peace_is_war(wp) ? "YES" : "NO");
    }
    war_peace_free(wp);

    /* 4. Euler buckling */
    printf("\n--- 4. Euler Buckling ---\n");
    EulerBuckling* eb = euler_buckling_create(1.0, 2.0, 0.5, 0);
    double Pcr = euler_buckling_critical_load(eb);
    printf("Beam: L=2.0 EI=0.5 Pcr=%.4f\n", Pcr);

    printf("Sweeping load from 0.5*Pcr to 2.0*Pcr:\n");
    for (int i = 0; i <= 10; i++) {
        eb->load = Pcr * (0.5 + 1.5 * i / 10.0);
        euler_buckling_simulate(eb, 2.0, 0.01);
        printf("  P/Pcr=%.2f  defl=%.4f  buckled=%s\n",
               eb->load / Pcr, eb->deflection,
               euler_buckling_has_buckled(eb) ? "YES" : "NO");
    }
    euler_buckling_free(eb);

    /* 5. Phase transition */
    printf("\n--- 5. van der Waals Phase Transition ---\n");
    /* CO2 parameters: a=3.64, b=0.0427 */
    PhaseTransition* pt = phase_transition_create(250.0, 0.1, 3.64, 0.0427);
    phase_transition_print(pt);

    printf("\nIsotherm sampling (T=250K, Tc=%.1fK):\n", pt->Tc);
    double V_vals[5], P_vals[5];
    phase_transition_isotherm(pt, V_vals, P_vals, 5);
    for (int i = 0; i < 5; i++) {
        printf("  V=%.4f  P=%.4f\n", V_vals[i], P_vals[i]);
    }

    /* 6. Ship stability */
    printf("\n--- 6. Ship Stability ---\n");
    ShipStability* ss = ship_stability_create(0.3, 12.0, 0.4);
    printf("Stable ship (GM=0.3): "); ship_stability_print(ss);
    ship_stability_free(ss);

    ShipStability* ss2 = ship_stability_create(-0.1, 12.0, 0.4);
    printf("Unstable ship (GM=-0.1): "); ship_stability_print(ss2);
    printf("  Capsize risk (30 deg): %s\n",
           ship_stability_capsize_check(ss2, 0.52) ? "HIGH" : "LOW");
    ship_stability_free(ss2);

    printf("\n=======================================\n");
    printf("  Demo complete.\n");
    printf("=======================================\n");
    return 0;
}