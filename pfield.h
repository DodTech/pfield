#ifndef PFIELDS_H
#define PFIELDS_H

#include<stdlib.h>
#include<math.h>

#include "pfield.h"
#include "kdtree.h"

typedef struct _pfield pfield;

pfield *pf_new(const unsigned short dim, const float a_coeff, 
    const float r_coeff, const float r_range);

void pf_free(pfield *pf);

int set_attractor(pfield *pf, const float q[]);

int new_repulsor(pfield *pf, const float q[]);

float pfield_energy(pfield *pf, const float q[]);

void pfield_force(pfield *pf, float g[], const float q[], const float step);

#endif // PFIELDS_H
