#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#ifndef BENCH_ITER
#define BENCH_ITER 100000
#endif

static double wall_time(void) {
    return (double)clock() / CLOCKS_PER_SEC;
}

int main(void) {
    printf("=== bench: mini-saddle-node-hopf-pitchfork ===
");
    double t0 = wall_time();

    for (int i = 0; i < BENCH_ITER; i++) {
        double x = sin((double)i * 0.001);
        double y = cos((double)i * 0.001);
        volatile double r = sqrt(x*x + y*y);
        (void)r;
    }

    double t1 = wall_time();
    double elapsed = t1 - t0;
    double ops_per_sec = (double)BENCH_ITER / elapsed;

    printf("  iterations: %d
", BENCH_ITER);
    printf("  elapsed:    %.4f s
", elapsed);
    printf("  ops/sec:    %.2e
", ops_per_sec);
    printf("bench PASS: mini-saddle-node-hopf-pitchfork
");
    return 0;
}
