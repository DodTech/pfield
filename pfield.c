#include<stdlib.h>
#include<math.h>

#include "pfield.h"
#include "kdtree.h"

struct _pfield {
    // types prefer space
    unsigned short dimensions;
    float a_coeff;
    float r_coeff;
    float r_range;
    float *attractor;  // array of coordinates
    void  *repulsors;  // kdtree storing coordinates
};

pfield *pf_new(const unsigned short dim, const float a_coeff, 
    const float r_coeff, const float r_range)
{
    pfield *pf = malloc(sizeof(pfield));
    pf->dimensions = dim;
    pf->a_coeff = a_coeff;
    pf->r_coeff = r_coeff;
    pf->r_range = r_range;
    pf->attractor = malloc(sizeof(float)*dim);
    pf->repulsors = kd_create(dim);

    return pf;
}

void pf_free(pfield *pf)
{
    free(pf->attractor);
    kd_free(pf->repulsors);
    free(pf);
}

int set_attractor(pfield *pf, const float q[])
{
    for(unsigned short i = 0; i < pf->dimensions; i++) {
        pf->attractor[i] = q[i];
    }
    
    return 0;
}

int new_repulsor(pfield *pf, const float q[])
{
    if(kd_insertf(pf->repulsors, q, NULL) != 0) {
        return -1;
    }
    
    return 0;
}

float dist_sq(const float a[], const float b[], unsigned short d)
{
    float dsq = 0;

    for(unsigned short i = 0; i < d; i++) {
        dsq += (a[i]-b[i])*(a[i]-b[i]);
    }

    return dsq;
}

float pfield_energy(pfield *pf, const float q[])
{
    struct kdres *kdr;
    float p[pf->dimensions];
    float potential;
    float dist;    
    float dsq;
    
    /* Attraction */
    dsq = dist_sq(pf->attractor, q, pf->dimensions);    
    potential = .5*pf->a_coeff*dsq;
    /* Repulsion */
    kdr = kd_nearest_rangef(pf->repulsors, q, pf->r_range);
    while(!kd_res_end(kdr)) {
        kd_res_itemf(kdr, p);
        dist = sqrt(dist_sq(p, q, pf->dimensions));
        potential += .5*pf->r_coeff*pow(1/dist-1/pf->r_range, 2);
        kd_res_next(kdr);
    }
    kd_res_free(kdr);

    return potential;
}

void pfield_force(pfield *pf, float g[], const float q[], const float step)
{
    float p[pf->dimensions];

    for(unsigned short i = 0; i < pf->dimensions; i++) {
        g[i] = 0;
        p[i] = q[i] + step;
        g[i] += pfield_energy(pf, p);
        p[i] = q[i] - step;
        g[i] -= pfield_energy(pf, p);
        p[i] = q[i];
        g[i] /= 2*step;
    }
}

