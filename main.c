#include <stdio.h>
#include "pfields.h"

int main(void)
{
    float o[] = {1.0, 0.0};
    float a[] = {4.0, 4.0};
    float c[] = {8.0, 8.0};
    float g[] = {1.0, 1.0};

    pfield *f = pf_new(2, .5, .5, 5);
    new_repulsor(f, a);
    set_attractor(f, c);
    
    printf("%.6f\n", pfield_energy(f, o));

    while(pfield_energy(f, o) > 1e-2) {
        pfield_force(f, g, o, .001);
        for(unsigned short j = 0; j < 2; j++) {
            o[j] -= g[j];
            printf("%9.6f ", o[j]);
        }
        printf("\n");
    }

    printf("%.6f\n", pfield_energy(f, o));

    pf_free(f);
    
    return 0;
}
