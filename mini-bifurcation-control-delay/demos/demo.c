#include <stdio.h>
#include "../include/bcd_core.h"

int main(void) {
    printf("=== Bifurcation Control Demo ===
");
    printf("Washout: dw/dt = x - w*w, u = -k*(x - w*w)
");
    printf("DFC (Pyragas): u(t) = K*[y(t-tau) - y(t)]
");
    printf("Amplitude death: coupling-induced quenching
");
    return 0;
}
