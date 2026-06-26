#include "bifurcation_detection.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

BifurcationScanner* bscan_create(void) {
    BifurcationScanner* bs = calloc(1, sizeof(BifurcationScanner));
    if (bs) {
        bs->sn = sn_create();
        bs->hb = hopf_create();
        bs->pb = pf_create();
    }
    return bs;
}

void bscan_free(BifurcationScanner* bs) {
    if (!bs) return;
    for (int i = 0; i < bs->n; i++) {
        free(bs->events[i].fp);
        eigen_free(bs->events[i].eigs);
    }
    free(bs->events);
    free(bs->param_scan);
    sn_free(bs->sn);
    hopf_free(bs->hb);
    pf_free(bs->pb);
    free(bs);
}

int bscan_add_event(BifurcationScanner* bs, BifurcationType t,
    double p, double* fp, int n) {
    if (!bs) return -1;
    if (bs->n >= bs->cap) {
        int nc = (bs->cap == 0) ? 16 : bs->cap * 2;
        BifurcationEvent* ne = realloc(bs->events,
            (size_t)nc * sizeof(BifurcationEvent));
        if (!ne) return -1;
        bs->events = ne;
        bs->cap = nc;
    }
    BifurcationEvent* ev = &bs->events[bs->n];
    ev->type = t;
    ev->param = p;
    ev->fp = malloc((size_t)n * sizeof(double));
    memcpy(ev->fp, fp, (size_t)n * sizeof(double));
    ev->n = n;
    ev->eigs = NULL;
    return bs->n++;
}

int bscan_detect_all(BifurcationScanner* bs, ODEFunc f,
    double* x0, double* params, int n, int pidx,
    double pmin, double pmax, int nsteps,
    int max_iter, double tol) {
    if (!bs || !f) return -1;
    bs->n_scan = nsteps;
    bs->param_scan = malloc((size_t)nsteps * sizeof(double));
    double dp = (pmax - pmin) / (double)(nsteps - 1);
    double x[8];
    memcpy(x, x0, (size_t)n * sizeof(double));

    for (int i = 0; i < nsteps; i++) {
        double p = pmin + (double)i * dp;
        bs->param_scan[i] = p;
        params[pidx] = p;

        int iter = fp_newton(f, x, params, n, max_iter, tol,
            &(FixedPoint){.point = x, .n = n});
        if (iter < 0) continue;

        EigenSpectrum* e = eigen_compute_jacobian(f, x, params, n, tol);
        if (!e) continue;

        if (sn_is_saddle_node(e)) {
            bscan_add_event(bs, BIF_SN, p, x, n);
        } else if (hopf_is_hopf(e, tol)) {
            bscan_add_event(bs, BIF_HOPF, p, x, n);
        } else if (pf_is_pitchfork(e, f, n)) {
            bscan_add_event(bs, BIF_PITCHFORK, p, x, n);
        }
        eigen_free(e);
    }
    return bs->n;
}

int bscan_count_type(const BifurcationScanner* bs, BifurcationType t) {
    if (!bs) return 0;
    int c = 0;
    for (int i = 0; i < bs->n; i++)
        if (bs->events[i].type == t) c++;
    return c;
}

BifurcationEvent* bscan_get_event(const BifurcationScanner* bs, int idx) {
    if (!bs || idx < 0 || idx >= bs->n) return NULL;
    return &bs->events[idx];
}

const char* bscan_type_name(BifurcationType t) {
    switch (t) {
        case BIF_SN: return "Saddle-Node";
        case BIF_HOPF: return "Hopf";
        case BIF_PITCHFORK: return "Pitchfork";
        case BIF_TRANS: return "Transcritical";
        default: return "None";
    }
}

void bscan_print(const BifurcationScanner* bs) {
    if (!bs) { printf("BifurcationScanner: NULL\n"); return; }
    printf("=== Bifurcation Scanner: %d events ===\n", bs->n);
    printf("  SN: %d  Hopf: %d  PF: %d  Trans: %d\n",
        bscan_count_type(bs, BIF_SN),
        bscan_count_type(bs, BIF_HOPF),
        bscan_count_type(bs, BIF_PITCHFORK),
        bscan_count_type(bs, BIF_TRANS));
    for (int i = 0; i < bs->n && i < 15; i++)
        printf("  [%d] p=%.8f  %s\n",
            i, bs->events[i].param,
            bscan_type_name(bs->events[i].type));
}

/* Extended scanner utilities */
int bscan_export_csv(const BifurcationScanner* bs, const char* fn) {
    if (!bs || !fn) return -1;
    FILE* f = fopen(fn, "w");
    if (!f) return -1;
    fprintf(f, "index,parameter,type,x\n");
    for (int i = 0; i < bs->n; i++)
        fprintf(f, "%d,%.8f,%s,%.8f\n",
            i, bs->events[i].param,
            bscan_type_name(bs->events[i].type),
            bs->events[i].fp[0]);
    fclose(f);
    return 0;
}

double bscan_first_bifurcation_param(const BifurcationScanner* bs) {
    if (!bs || bs->n == 0) return 0.0;
    return bs->events[0].param;
}

int bscan_largest_event_index(const BifurcationScanner* bs) {
    if (!bs || bs->n == 0) return -1;
    return bs->n - 1;
}
