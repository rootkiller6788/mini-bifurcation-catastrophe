#include "zeeman_social.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Prison Riot Model (Zeeman 1976)
 *
 * Collective behavior of prisoners modeled as a cusp catastrophe.
 *
 * Control parameters:
 *   - Tension (a): accumulated grievance, alienation
 *   - Trigger (b): immediate provoking event
 *
 * State variable x: degree of riot/insurrection
 *   x < 0: calm, ordered behavior
 *   x > 0: riot, collective violence
 *
 * The cusp structure explains:
 *   - Why identical prisons with similar tension levels can have
 *     vastly different outcomes (bimodality)
 *   - Why a small triggering event can cause a sudden riot
 *     (catastrophic jump)
 *   - Why once a riot starts, tension must drop significantly
 *     below the onset level to restore order (hysteresis)
 *
 * Reference: Zeeman (1976), "Prison Riots", in Catastrophe Theory:
 * Selected Papers 1972-1977.
 * ============================================================== */

PrisonRiotModel* prison_riot_create(double tension, double trigger) {
    PrisonRiotModel* p = calloc(1, sizeof(PrisonRiotModel));
    if (!p) return NULL;
    p->state = -0.5;    /* start calm */
    p->tension = tension;
    p->trigger = trigger;
    /* Map to cusp parameters:
     * a_factor = -(tension - tension_critical) -> negative when tense
     * b_factor = trigger - trigger_baseline */
    p->a_factor = -(tension - 0.5);   /* 0.5 is neutral tension */
    p->b_factor = trigger - 0.3;      /* 0.3 is baseline trigger */
    p->t = 0.0;
    return p;
}

void prison_riot_free(PrisonRiotModel* p) {
    free(p);
}

double prison_riot_potential(const PrisonRiotModel* p, double x) {
    return cusp_potential(x, p->a_factor, p->b_factor);
}

CuspEquilibrium* prison_riot_equilibrium(const PrisonRiotModel* p) {
    return cusp_equilibrium(p->a_factor, p->b_factor);
}

void prison_riot_step(PrisonRiotModel* p, double dt) {
    /* Update cusp parameters from tension/trigger */
    p->a_factor = -(p->tension - 0.5);
    p->b_factor = p->trigger - 0.3;

    /* Gradient dynamics with social inertia */
    double force = cusp_force(p->state, p->a_factor, p->b_factor);
    double noise = 0.0;
    /* Add small stochastic component for more realistic behavior */
    noise = 0.02 * ((double)rand() / RAND_MAX - 0.5);
    p->state += (force + noise) * dt;
    p->t += dt;

    if (p->state > 2.0)  p->state = 2.0;
    if (p->state < -2.0) p->state = -2.0;

    /* Gradual tension escalation when in riot state */
    if (p->state > 0.5) {
        p->tension += 0.01 * dt;  /* riot increases tension */
    } else if (p->state < -0.5) {
        p->tension -= 0.02 * dt;  /* calm decreases tension */
        if (p->tension < 0.0) p->tension = 0.0;
    }
    if (p->tension > 1.0) p->tension = 1.0;
}

void prison_riot_simulate(PrisonRiotModel* p, double duration, double dt) {
    int steps = (int)(duration / dt);
    for (int i = 0; i < steps; i++) {
        prison_riot_step(p, dt);
    }
}

double prison_riot_critical_tension(const PrisonRiotModel* p) {
    /* Critical tension where riot becomes possible (entering cusp region)
     * Need a_factor < 0, so -(tension - 0.5) < 0 => tension > 0.5
     * And inside cusp: 4*a^3 + 27*b^2 < 0
     * a = -(T - 0.5), b = trigger - 0.3
     * 4*(-(T-0.5))^3 + 27*(trigger-0.3)^2 < 0
     * -4*(T-0.5)^3 + 27*(trigger-0.3)^2 < 0
     * 4*(T-0.5)^3 > 27*(trigger-0.3)^2
     * T > 0.5 + cbrt(27*(trigger-0.3)^2 / 4) */
    double b_sq = (p->trigger - 0.3) * (p->trigger - 0.3);
    return 0.5 + pow(27.0 * b_sq / 4.0, 1.0/3.0);
}

bool prison_riot_is_riot(const PrisonRiotModel* p) {
    return p->state > 0.3;
}

void prison_riot_hysteresis(const PrisonRiotModel* p, double t_min, double t_max,
                             int steps, HysteresisLoop* loop) {
    /* Vary tension while keeping trigger fixed */
    PrisonRiotModel* pw = prison_riot_create(t_min, p->trigger);
    if (!pw) return;

    /* Forward sweep: increasing tension */
    for (int i = 0; i <= steps; i++) {
        pw->tension = t_min + (t_max - t_min) * i / steps;
        pw->a_factor = -(pw->tension - 0.5);
        pw->b_factor = pw->trigger - 0.3;
        /* Find equilibrium */
        CuspEquilibrium* eq = cusp_equilibrium(pw->a_factor, pw->b_factor);
        if (eq && eq->n_roots > 0) {
            pw->state = eq->x[0];
        }
        cusp_equilibrium_free(eq);
        if (loop) {
            loop->a_fwd[i] = pw->tension;
            loop->x_fwd[i] = pw->state;
        }
    }

    /* Backward sweep: decreasing tension */
    for (int i = 0; i <= steps; i++) {
        pw->tension = t_max + (t_min - t_max) * i / steps;
        pw->a_factor = -(pw->tension - 0.5);
        pw->b_factor = pw->trigger - 0.3;
        CuspEquilibrium* eq = cusp_equilibrium(pw->a_factor, pw->b_factor);
        if (eq && eq->n_roots > 0) {
            /* Stay on current sheet if possible */
            int best = 0;
            double min_dist = fabs(eq->x[0] - pw->state);
            for (int j = 1; j < eq->n_roots; j++) {
                double dist = fabs(eq->x[j] - pw->state);
                if (dist < min_dist) { min_dist = dist; best = j; }
            }
            pw->state = eq->x[best];
        }
        cusp_equilibrium_free(eq);
        if (loop) {
            loop->a_bwd[i] = pw->tension;
            loop->x_bwd[i] = pw->state;
        }
    }
    prison_riot_free(pw);
}

void prison_riot_print(const PrisonRiotModel* p) {
    printf("PrisonRiot: state=%.4f tension=%.3f trigger=%.3f a=%.3f b=%.3f t=%.2f\n",
           p->state, p->tension, p->trigger, p->a_factor, p->b_factor, p->t);
    printf("  Status: %s  CriticalTension: %.3f  Bistable: %s\n",
           prison_riot_is_riot(p) ? "RIOT" : "CALM",
           prison_riot_critical_tension(p),
           cusp_is_bistable(p->a_factor, p->b_factor) ? "YES" : "NO");
}

/* ==============================================================
 * Stock Market Model (Zeeman 1974)
 *
 * Financial market behavior modeled as a cusp catastrophe.
 *
 * Control parameters:
 *   - Speculation (a): degree of speculative excess, momentum trading
 *   - Fundamental valuation (b): deviation from fundamental value
 *
 * State variable x: market regime
 *   x < 0: bear market
 *   x > 0: bull market
 *
 * The cusp explains:
 *   - Sudden market crashes without commensurate news (catastrophe)
 *   - Why markets can stay irrationally high then crash (hysteresis)
 *   - Divergence between fundamentally similar stocks
 *
 * Reference: Zeeman (1974), "On the Unstable Behavior of
 * Stock Exchanges", Journal of Mathematical Economics.
 * ============================================================== */

StockMarketModel* stock_market_create(double price, double speculation,
                                       double fundamental, double volatility) {
    StockMarketModel* s = calloc(1, sizeof(StockMarketModel));
    if (!s) return NULL;
    s->price = price;
    s->speculation = speculation;
    s->fundamental = fundamental;
    s->volatility = volatility;
    /* Map to cusp parameters */
    s->a_factor = -(speculation - 0.5);
    s->b_factor = (price - fundamental) / fundamental;  /* relative mispricing */
    s->t = 0.0;
    return s;
}

void stock_market_free(StockMarketModel* s) {
    free(s);
}

double stock_market_potential(const StockMarketModel* s, double x) {
    return cusp_potential(x, s->a_factor, s->b_factor);
}

CuspEquilibrium* stock_market_equilibrium(const StockMarketModel* s) {
    return cusp_equilibrium(s->a_factor, s->b_factor);
}

void stock_market_step(StockMarketModel* s, double dt) {
    /* Update cusp parameters */
    s->a_factor = -(s->speculation - 0.5);
    s->b_factor = (s->price - s->fundamental) / fmax(s->fundamental, 1.0);

    /* Price dynamics with cusp-driven regime switching */
    double force = cusp_force(s->price / s->fundamental - 1.0,
                               s->a_factor, s->b_factor);
    double noise = s->volatility * ((double)rand() / RAND_MAX - 0.5);
    s->price += (force * 10.0 + noise) * s->fundamental * dt;
    s->t += dt;

    if (s->price < 0.1 * s->fundamental) s->price = 0.1 * s->fundamental;

    /* Speculation decays slowly, but can spike */
    s->speculation += (-0.05 * s->speculation
                       + 0.01 * (s->price - s->fundamental) / s->fundamental) * dt;
    if (s->speculation > 1.0) s->speculation = 1.0;
    if (s->speculation < 0.0) s->speculation = 0.0;
}

void stock_market_simulate(StockMarketModel* s, double duration, double dt) {
    int steps = (int)(duration / dt);
    for (int i = 0; i < steps; i++) {
        stock_market_step(s, dt);
    }
}

bool stock_market_crash_detect(const StockMarketModel* s, double threshold) {
    /* Crash: price drops below fundamental by more than threshold AND
     * the system is in the cusp region (potential for sudden jump) */
    double deviation = (s->price - s->fundamental) / fmax(s->fundamental, 1.0);
    return deviation < -threshold && cusp_is_bistable(s->a_factor, s->b_factor);
}

double stock_market_crash_magnitude(const StockMarketModel* s) {
    /* Magnitude of potential crash = distance between stable equilibria */
    CuspEquilibrium* eq = cusp_equilibrium(s->a_factor, s->b_factor);
    double mag = 0.0;
    if (eq && eq->n_stable >= 2) {
        /* Find two stable equilibria */
        int s1 = -1, s2 = -1;
        for (int i = 0; i < eq->n_roots; i++) {
            if (eq->stability[i] == STABILITY_STABLE) {
                if (s1 < 0) s1 = i;
                else { s2 = i; break; }
            }
        }
        if (s2 >= 0) {
            mag = fabs(eq->x[s1] - eq->x[s2]);
        }
    }
    cusp_equilibrium_free(eq);
    return mag * s->fundamental;
}

void stock_market_hysteresis(const StockMarketModel* s, double spec_min,
                              double spec_max, int steps, HysteresisLoop* loop) {
    for (int i = 0; i <= steps; i++) {
        double spec = spec_min + (spec_max - spec_min) * i / steps;
        double a = -(spec - 0.5);
        double b = (s->price - s->fundamental) / fmax(s->fundamental, 1.0);
        CuspEquilibrium* eq = cusp_equilibrium(a, b);
        double x = (eq && eq->n_roots > 0) ? eq->x[0] : 0.0;
        cusp_equilibrium_free(eq);
        if (loop) {
            loop->a_fwd[i] = spec;
            loop->x_fwd[i] = x;
        }
    }
    for (int i = 0; i <= steps; i++) {
        double spec = spec_max + (spec_min - spec_max) * i / steps;
        double a = -(spec - 0.5);
        double b = (s->price - s->fundamental) / fmax(s->fundamental, 1.0);
        CuspEquilibrium* eq = cusp_equilibrium(a, b);
        double x = (eq && eq->n_roots > 0) ? eq->x[0] : 0.0;
        /* For backward: stay on upper sheet if possible */
        if (eq && eq->n_roots > 1) x = eq->x[eq->n_roots - 1];
        cusp_equilibrium_free(eq);
        if (loop) {
            loop->a_bwd[i] = spec;
            loop->x_bwd[i] = x;
        }
    }
}

void stock_market_print(const StockMarketModel* s) {
    printf("StockMarket: price=%.2f spec=%.3f fundamental=%.2f vol=%.3f t=%.2f\n",
           s->price, s->speculation, s->fundamental, s->volatility, s->t);
    printf("  Dev: %.1f%%  CrashRisk: %s  CrashMag: %.2f\n",
           100.0 * (s->price - s->fundamental) / s->fundamental,
           stock_market_crash_detect(s, 0.2) ? "HIGH" : "LOW",
           stock_market_crash_magnitude(s));
}

/* ==============================================================
 * War/Peace Model (Zeeman 1979)
 *
 * International conflict modeled as a cusp catastrophe.
 *
 * Control parameters:
 *   - Threat perception (a): how threatening the adversary appears
 *   - Cost-benefit ratio (b): net benefit of going to war
 *
 * State variable x: diplomatic/military state
 *   x < 0: peace / diplomacy
 *   x > 0: war / armed conflict
 *
 * The cusp explains:
 *   - Sudden outbreak of war (catastrophic jump from peace to war)
 *   - Why peace can persist even under high threat (hysteresis)
 *   - Why de-escalation is harder than escalation (asymmetry)
 *
 * Reference: Zeeman (1979), "A Model for Conflicts", in
 * Catastrophe Theory: Selected Papers.
 * ============================================================== */

WarPeaceModel* war_peace_create(double threat, double cost) {
    WarPeaceModel* w = calloc(1, sizeof(WarPeaceModel));
    if (!w) return NULL;
    w->state = -1.0;    /* start at peace */
    w->threat = threat;
    w->cost = cost;
    w->a_factor = -(threat - 0.5);   /* high threat -> negative a -> bistable */
    w->b_factor = 0.5 - cost;        /* low cost -> positive b -> war bias */
    w->t = 0.0;
    return w;
}

void war_peace_free(WarPeaceModel* w) {
    free(w);
}

double war_peace_potential(const WarPeaceModel* w, double x) {
    return cusp_potential(x, w->a_factor, w->b_factor);
}

CuspEquilibrium* war_peace_equilibrium(const WarPeaceModel* w) {
    return cusp_equilibrium(w->a_factor, w->b_factor);
}

void war_peace_step(WarPeaceModel* w, double dt) {
    w->a_factor = -(w->threat - 0.5);
    w->b_factor = 0.5 - w->cost;

    /* Diplomatic/military dynamics */
    double force = cusp_force(w->state, w->a_factor, w->b_factor);
    w->state += force * dt;
    w->t += dt;

    if (w->state > 2.0)  w->state = 2.0;
    if (w->state < -2.0) w->state = -2.0;

    /* Feedback: war state increases threat perception */
    if (w->state > 0.0) {
        w->threat += 0.01 * dt;
        if (w->threat > 1.0) w->threat = 1.0;
    }
}

void war_peace_simulate(WarPeaceModel* w, double duration, double dt) {
    int steps = (int)(duration / dt);
    for (int i = 0; i < steps; i++) {
        war_peace_step(w, dt);
    }
}

double war_peace_threshold(const WarPeaceModel* w) {
    /* Threshold threat level where war becomes possible */
    if (w->a_factor >= 0.0) return 1.0; /* no bistability */
    return 0.5 + pow(27.0 * w->b_factor * w->b_factor / 4.0, 1.0/3.0);
}

bool war_peace_is_war(const WarPeaceModel* w) {
    return w->state > 0.0;
}

void war_peace_hysteresis(const WarPeaceModel* w, double threat_min,
                           double threat_max, int steps, HysteresisLoop* loop) {
    for (int i = 0; i <= steps; i++) {
        double threat = threat_min + (threat_max - threat_min) * i / steps;
        double a = -(threat - 0.5);
        double b = w->b_factor;
        CuspEquilibrium* eq = cusp_equilibrium(a, b);
        double x = (eq && eq->n_roots > 0) ? eq->x[0] : 0.0;
        cusp_equilibrium_free(eq);
        if (loop) { loop->a_fwd[i] = threat; loop->x_fwd[i] = x; }
    }
    for (int i = 0; i <= steps; i++) {
        double threat = threat_max + (threat_min - threat_max) * i / steps;
        double a = -(threat - 0.5);
        double b = w->b_factor;
        CuspEquilibrium* eq = cusp_equilibrium(a, b);
        double x;
        if (eq && eq->n_roots > 1) x = eq->x[eq->n_roots - 1]; /* upper sheet */
        else if (eq && eq->n_roots > 0) x = eq->x[0];
        else x = 0.0;
        cusp_equilibrium_free(eq);
        if (loop) { loop->a_bwd[i] = threat; loop->x_bwd[i] = x; }
    }
}

void war_peace_print(const WarPeaceModel* w) {
    printf("WarPeace: state=%.4f threat=%.3f cost=%.3f a=%.3f b=%.3f t=%.2f\n",
           w->state, w->threat, w->cost, w->a_factor, w->b_factor, w->t);
    printf("  Status: %s  Threshold: %.3f  Bistable: %s\n",
           war_peace_is_war(w) ? "WAR" : "PEACE",
           war_peace_threshold(w),
           cusp_is_bistable(w->a_factor, w->b_factor) ? "YES" : "NO");
}