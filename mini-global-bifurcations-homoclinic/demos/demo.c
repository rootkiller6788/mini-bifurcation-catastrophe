#include <stdio.h>
#include "../include/homoclinic_core.h"

int main(void) {
    printf("=== Homoclinic Bifurcation Demo ===
");
    printf("Homoclinic: trajectory returns to saddle as t->+/-inf
");
    printf("Melnikov: M(t0) = integral of wedge product
");
    printf("Shilnikov: saddle-focus homoclinic -> chaos
");
    return 0;
}
