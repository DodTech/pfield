#include<stdlib.h>
#include<math.h>

#include "pfield.h"
#include "kdtree.h"

/*
 *      struct: _pfield 
 *      ---------------
 *      Artificial potential field. Only one attractor can be present. Many 
 *      repulsors are stored in a kd-tree. Types generally prefer space.
 *      
 *      dim: number of dimensions
 *      a_coeff: attractive coefficient (eg 2.0)
 *      r_coeff: repulsive coefficient (eg 1.5)
 *      a_range: threshold of quadratic attraction (eg 5)
 *      r_range: threshold of repulsion (eg .5)
 */

struct _pfield {
    unsigned short dimensions;
    float a_coeff;
    float r_coeff;
    float a_range;
    float r_range;
    float *attractor;
    void  *repulsors;
};

/*
 *      function: pf_new
 *      ---------------- 
 *      Allocates and returns new potential field
 *      
 *      dim: number of dimensions
 *      a_coeff: attractive coefficient (eg 2.0)
 *      r_coeff: repulsive coefficient (eg 1.5)
 *      a_range: threshold of quadratic attraction (eg 5)
 *      r_range: threshold of repulsion (eg .5)
 *      
 *      returns: new potential field
 */

pfield *pf_new(const unsigned short dim, 
    const float a_coeff, const float r_coeff, 
    const float a_range, const float r_range)
{
    pfield *pf = malloc(sizeof(pfield));
    pf->dimensions = dim;
    pf->a_coeff = a_coeff;
    pf->r_coeff = r_coeff;
    pf->a_range = a_range;
    pf->r_range = r_range;
    pf->attractor = malloc(sizeof(float)*dim);
    pf->repulsors = kd_create(dim);

    return pf;
}

/*
 *      function: pf_free
 *      -----------------
 *      Frees kd-tree, array, and struct
 *      
 *      pf: potential field to be freed
 */

void pf_free(pfield *pf)
{
    free(pf->attractor);
    kd_free(pf->repulsors);
    free(pf);
}

/*
 *      function: set_attractor
 *      -----------------------
 *      Set attractor position
 *
 *      pf: potential field
 *      q: position of attractor
 *
 *      returns: 0 if successful
 */

int set_attractor(pfield *pf, const float q[])
{
    for(unsigned short i = 0; i < pf->dimensions; i++) {
        pf->attractor[i] = q[i];
    }
    
    return 0;
}

/*
 *      function: new_repulsor
 *      ----------------------
 *      Add new repulsor position to kd-tree
 *
 *      pf: potential field
 *      q: position of repulsor
 *
 *      returns: 0 if successful    
 */

int new_repulsor(pfield *pf, const float q[])
{
    if(kd_insertf(pf->repulsors, q, NULL) != 0) {
        return -1;
    }
    
    return 0;
}

/*
 *      function: distance
 *      ------------------
 *      Calculates euclidean distance between two coordinates
 *
 *      a: coordinate
 *      b: coordinate
 *      d: dimensions of the coordinates
 *
 *      returns: distance
 */

float distance(const float a[], const float b[], unsigned short d)
{
    float sq_sum = 0;
    float dist;

    for(unsigned short i = 0; i < d; i++) {
        sq_sum += pow(a[i]-b[i], 2);
    }

    dist = sqrt(sq_sum);

    return dist;
}

/*
 *      function: attraction_energy
 *      ---------------------------
 *      Calculates attraction energy at point
 *
 *      pf: potential field
 *      q: current position
 *
 *      returns: attraction energy
 */

float attraction_energy(pfield *pf, const float q[])
{
    float dist, attraction;

    dist = distance(pf->attractor, q, pf->dimensions);
    if(dist < pf->a_range) {    
        attraction = .5*pf->a_coeff*pow(dist, 2);
    } else {
        attraction = pf->a_range*pf->a_coeff*(dist-.5*pf->a_range);
    }

    return attraction;
}

/*
 *      function: repulsion_energy
 *      ---------------------------
 *      Calculates repulsion energy at point
 *
 *      pf: potential field
 *      q: current position
 *
 *      returns: attraction energy
 */

float repulsion_energy(pfield *pf, const float q[])
{
    struct kdres *kdr;
    float p[pf->dimensions];
    float repulsion;
    float dist;

    repulsion = 0;
    kdr = kd_nearest_rangef(pf->repulsors, q, pf->r_range);
    while(!kd_res_end(kdr)) {
        kd_res_itemf(kdr, p);
        dist = distance(p, q, pf->dimensions);
        repulsion += .5*pf->r_coeff*pow(1/dist-1/pf->r_range, 2);
        kd_res_next(kdr);
    }
    kd_res_free(kdr);

    return repulsion;
}

/*
 *      function: pfield_energy
 *      -----------------------
 *      Calculates the potential energy at the given position. Energy is 
 *      the sum of attraction and repulsion
 *
 *      pf: potential field
 *      q: current position
 *
 *      returns: energy at position
 */

float pfield_energy(pfield *pf, const float q[])
{
    float potential;
    
    potential = attraction_energy(pf, q) + repulsion_energy(pf, q);

    return potential;
}

/*
 *      function: pfield_force
 *      ----------------------
 *      Calculates the force at the given position. Found with 3 point
 *      midpoint numerical gradient
 *
 *      pf: potential field
 *      g: gradient, value is changed
 *      q: current position
 *      step: step size for gradient
 */

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

/*
 *      function: pfield_move
 *      ---------------------
 *      Moves along the normalize gradient of the potential field
 *
 *      pf: potential field
 *      q: current position, value is changed
 *      step: step size for move
 *
 *      returns: 0 if successful
 */

int pfield_move(pfield *pf, float q[], const float step)
{
    float g[pf->dimensions];
    float ssum, norm;

    pfield_force(pf, g, q, step);

    ssum = 0;
    for(unsigned short i = 0; i < pf->dimensions; i++) {
        ssum += pow(g[i], 2);
    }
    norm = sqrt(ssum);

    if(norm == 0) {
        return -1;
    }

    for(unsigned short i = 0; i < pf->dimensions; i++) {
        q[i] -= step*g[i]/norm;
    }

    return 0;
}
