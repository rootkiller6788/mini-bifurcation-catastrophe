#ifndef ZEEMAN_BIOLOGY_H
#define ZEEMAN_BIOLOGY_H
#include "zeeman_core.h"

/* Zeeman's Biological Applications of Catastrophe Theory */

typedef struct {
    double x;       /* muscle fiber length */
    double y;       /* electrochemical activity */
    double a0;      /* baseline normal factor */
    double b0;      /* baseline splitting factor */
    double epsilon; /* time-scale separation */
    double I;       /* pacemaker stimulus */
    double t;
} HeartbeatModel;

typedef struct {
    double v;       /* membrane potential */
    double w;       /* recovery variable */
    double a_param; /* activation threshold */
    double b_param; /* recovery rate */
    double I_ext;   /* external stimulus */
    double t;
} NerveImpulse;

typedef struct {
    double x;           /* cell state */
    double y;           /* gene regulatory activity */
    double morphogen;   /* external morphogen concentration */
    double threshold;   /* critical level */
    double hill_coef;   /* cooperativity */
    double t;
} CellDiffModel;

/* Heartbeat model API (Zeeman 1972, 1977) */
HeartbeatModel* heartbeat_create(double a0, double b0, double epsilon, double I);
void            heartbeat_free(HeartbeatModel* h);
void            heartbeat_rhs(const HeartbeatModel* h, double* dx, double* dy);
void            heartbeat_step(HeartbeatModel* h, double dt);
void            heartbeat_simulate(HeartbeatModel* h, double duration, double dt);
bool            heartbeat_is_systole(const HeartbeatModel* h);
double          heartbeat_period(HeartbeatModel* h, double duration, double dt);
void            heartbeat_find_fixed_points(const HeartbeatModel* h, CuspEquilibrium* eq_out);
void            heartbeat_print(const HeartbeatModel* h);

/* Nerve impulse API (Zeeman 1973) */
NerveImpulse* nerve_impulse_create(double a_param, double b_param, double I_ext);
void          nerve_impulse_free(NerveImpulse* n);
void          nerve_impulse_rhs(const NerveImpulse* n, double* dv, double* dw);
void          nerve_impulse_step(NerveImpulse* n, double dt);
void          nerve_impulse_simulate(NerveImpulse* n, double duration, double dt);
double        nerve_impulse_threshold(const NerveImpulse* n);
bool          nerve_impulse_is_firing(const NerveImpulse* n);
int           nerve_impulse_spike_count(NerveImpulse* n, double duration, double dt);
void          nerve_impulse_print(const NerveImpulse* n);

/* Cell differentiation API (Thom 1975, Zeeman 1974) */
CellDiffModel* cell_diff_create(double morphogen, double threshold, double hill_coef);
void           cell_diff_free(CellDiffModel* c);
double         cell_diff_potential(const CellDiffModel* c, double x);
void           cell_diff_rhs(const CellDiffModel* c, double* dx, double* dy);
void           cell_diff_step(CellDiffModel* c, double dt);
void           cell_diff_simulate(CellDiffModel* c, double duration, double dt);
double         cell_diff_critical_morphogen(const CellDiffModel* c);
int            cell_diff_stable_states(const CellDiffModel* c);
void           cell_diff_print(const CellDiffModel* c);

#endif /* ZEEMAN_BIOLOGY_H */