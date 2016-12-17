#ifndef PFIELDS_H
#define PFIELDS_H

#include<stdlib.h>
#include<math.h>

#include "pfield.h"
#include "kdtree.h"

/*
 *      type: pfield 
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

typedef struct _pfield pfield;

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
    const float a_range, const float r_range);

/*
 *      function: pf_free
 *      -----------------
 *      Frees kd-tree, array, and struct
 *      
 *      pf: potential field to be freed
 */

void pf_free(pfield *pf);

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

int set_attractor(pfield *pf, const float q[]);

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

int new_repulsor(pfield *pf, const float q[]);

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

float pfield_energy(pfield *pf, const float q[]);

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

void pfield_force(pfield *pf, float g[], const float q[], const float step);

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

int pfield_move(pfield *pf, float q[], const float step);

#endif // PFIELDS_H
