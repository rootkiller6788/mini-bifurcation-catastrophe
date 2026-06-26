#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "homoclinic_core.h"
#include "homoclinic_orbit.h"
#include "homoclinic_shilnikov.h"
#include "homoclinic_melnikov.h"
#include "homoclinic_tangle.h"
#include "homoclinic_computation.h"
#define EPS 1e-9
int main(void){
    printf("=== Homoclinic Bifurcation Tests ===\n");
    HOMSystem* s=hom_create(3,4,0.01);assert(s);assert(s->dim==3);
    double x[]={1,1,1};hom_set_state(s,x);
    assert(fabs(hom_get_state(s,0)-1)<EPS);assert(fabs(hom_get_state(s,1)-1)<EPS);
    hom_set_param(s,0,10.0);assert(fabs(s->params[0]-10.0)<EPS);
    HOMSystem* sh=hom_create_shilnikov(0.1,1.0,-0.5,0.01,0.01);assert(sh);
    hom_step_rk4(sh);assert(sh->t>0);
    HOMOrbit* o=hom_orbit_create(100);assert(o);
    HOMConnection* c=hom_connection_create();assert(c);
    HOMShilnikovAnalysis* sa=hom_shilnikov_create(10);assert(sa);
    HOMTangle* ta=hom_tangle_create();assert(ta);
    HOMContinuation* co=hom_continuation_create(0,0,1,10);assert(co);
    HOMState g;g.dim=3;g.x[0]=0;g.x[1]=0;g.x[2]=0;
    HOMShooting* sho=hom_shooting_create(&g,3);assert(sho);
    hom_shooting_solve(sh,NULL,sho,1e-3,10);assert(hom_shooting_converged(sho));
    //assert(hom_bifurcation_codimension(sh,NULL)>=0);
    assert(!hom_is_global_bifurcation(sh,NULL,0.1));
    hom_free(s);hom_free(sh);hom_orbit_free(o);hom_connection_free(c);
    hom_shilnikov_free(sa);hom_tangle_free(ta);
    hom_continuation_free(co);hom_shooting_free(sho);
    printf("All tests passed.\n");return 0;
}
    assert(s->dim == 3); /* verify dimension */
