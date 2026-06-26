#ifndef ZEEMAN_SOCIAL_H
#define ZEEMAN_SOCIAL_H
#include "zeeman_core.h"

/* Zeeman's Social Science Applications of Catastrophe Theory */

typedef struct {
    double state;       /* collective behavior: calm <-> riot */
    double tension;     /* alienation, grievance level */
    double trigger;     /* immediate provoking factor */
    double a_factor;    /* derived normal factor */
    double b_factor;    /* derived splitting factor */
    double t;
} PrisonRiotModel;

typedef struct {
    double price;           /* current market price */
    double speculation;     /* speculative excess */
    double fundamental;     /* fundamental value anchor */
    double volatility;      /* noise amplitude */
    double a_factor;
    double b_factor;
    double t;
} StockMarketModel;

typedef struct {
    double state;       /* diplomatic state: peace <-> war */
    double threat;      /* perceived threat level */
    double cost;        /* estimated cost of war */
    double a_factor;
    double b_factor;
    double t;
} WarPeaceModel;

/* Prison riot API (Zeeman 1976) */
PrisonRiotModel* prison_riot_create(double tension, double trigger);
void             prison_riot_free(PrisonRiotModel* p);
double           prison_riot_potential(const PrisonRiotModel* p, double x);
CuspEquilibrium* prison_riot_equilibrium(const PrisonRiotModel* p);
void             prison_riot_step(PrisonRiotModel* p, double dt);
void             prison_riot_simulate(PrisonRiotModel* p, double duration, double dt);
double           prison_riot_critical_tension(const PrisonRiotModel* p);
bool             prison_riot_is_riot(const PrisonRiotModel* p);
void             prison_riot_hysteresis(const PrisonRiotModel* p, double t_min, double t_max, int steps, HysteresisLoop* loop);
void             prison_riot_print(const PrisonRiotModel* p);

/* Stock market API (Zeeman 1974) */
StockMarketModel* stock_market_create(double price, double speculation, double fundamental, double volatility);
void              stock_market_free(StockMarketModel* s);
double            stock_market_potential(const StockMarketModel* s, double x);
CuspEquilibrium*  stock_market_equilibrium(const StockMarketModel* s);
void              stock_market_step(StockMarketModel* s, double dt);
void              stock_market_simulate(StockMarketModel* s, double duration, double dt);
bool              stock_market_crash_detect(const StockMarketModel* s, double threshold);
double            stock_market_crash_magnitude(const StockMarketModel* s);
void              stock_market_hysteresis(const StockMarketModel* s, double spec_min, double spec_max, int steps, HysteresisLoop* loop);
void              stock_market_print(const StockMarketModel* s);

/* War/peace API (Zeeman 1979) */
WarPeaceModel* war_peace_create(double threat, double cost);
void           war_peace_free(WarPeaceModel* w);
double         war_peace_potential(const WarPeaceModel* w, double x);
CuspEquilibrium* war_peace_equilibrium(const WarPeaceModel* w);
void           war_peace_step(WarPeaceModel* w, double dt);
void           war_peace_simulate(WarPeaceModel* w, double duration, double dt);
double         war_peace_threshold(const WarPeaceModel* w);
bool           war_peace_is_war(const WarPeaceModel* w);
void           war_peace_hysteresis(const WarPeaceModel* w, double threat_min, double threat_max, int steps, HysteresisLoop* loop);
void           war_peace_print(const WarPeaceModel* w);

#endif /* ZEEMAN_SOCIAL_H */