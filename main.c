#include <stdio.h>
#include <stdlib.h>
#include "pfield.h"

int grid_energy(char name[], pfield *f, float a, float b, float step)
{
    float p[] = {0, 0};
    float e;

    FILE *fp = fopen(name, "w+");
    fprintf(fp, "# X\tY\tZ\n");

    for(float i = 0; i <= a; i+=step) {
        p[0] = i;
        for(float j = 0; j <= b; j+=step) {
            p[1] = j;
            e = pfield_energy(f, p);
            fprintf(fp, "%f\t%f\t%f\n", i, j, e);
        }
        fprintf(fp, "\n");
        printf("\rPercent Complete: %5.2f", 100*i/a);
    }
    printf("\rPercent Complete: 100.00\n");

    return fclose(fp);
}

int grid_force(char name[], pfield *f, float a, float b, float step)
{
    float p[] = {0, 0};
    float g[] = {0, 0};

    FILE *fp = fopen(name, "w+");
    fprintf(fp, "# X\tY\tZ1\tZ2\n");

    for(float i = 0; i <= a; i+=step) {
        p[0] = i;
        for(float j = 0; j <= b; j+=step) {
            p[1] = j;
            pfield_force(f, g, p, step);
            fprintf(fp, "%f\t%f\t%f\t%f\n", i, j, g[0], g[1]);
        }
        fprintf(fp, "\n");
        printf("\rPercent Complete: %.2f", 100*i / a);
    }
    printf("\n");

    return fclose(fp);
}

int path_traverse(char name[], pfield *f, int dim, 
    float p[], float tol, float step, float alpha)
{
    float e;
    float z = 0;

    FILE *fp = fopen(name, "w+");
    fprintf(fp, "# X\tY\tZ\n");
 
    do {
        if(pfield_move(f, p, step) != 0) {
            printf("Underflow in gradient norm. Stopped.");
            break;
        }
        for(unsigned short j = 0; j < dim; j++) {
            fprintf(fp, "%9.6f\t", p[j]);
        }
        
        e = pfield_energy(f, p);
        z = alpha*e + (1 - alpha)*z;

        if(fabs(z - e) < tol) {
            printf("\nFiltered and nominal converged. Stopped.");
            break;
        }

        printf("\rCurrent Energy: %2.6f %2.6f", z, e);
        fprintf(fp, "%9.6f\n", e);
    } while(e > tol);

    printf("\n");

    return fclose(fp);
}

int main(void)
{
    srand(1234);

    float g[] = {0.0, 0.0};
    float t[] = {5.0, 5.0};

    pfield *f = pf_new(2, 1.5, 2.0, 5.0, 0.5);

    set_attractor(f, t);
    for(int i = 0; i < 10; i++) {
        t[0] = (rand() % 4000)/1000.0 + 1;
        t[1] = (rand() % 4000)/1000.0 + 1;
        new_repulsor(f, t);
    }

    grid_energy("out/energy.tsv", f, 6, 6, .01);

    path_traverse("out/path.tsv", f, 2, g, .001, .0001, .1);

    pf_free(f);

    return 0;
}
