#include<stdio.h>
#include<math.h>
#include<libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"


/*---------------------------------------
  enuc computes nuclear repulsion energy
 ---------------------------------------*/

void compute_enuc()
{
  int i, j;
  double r = 0.0;
  double E = 0.0;

  if(Molecule.num_atoms > 1)
    for(i=1; i<Molecule.num_atoms; i++)
      for(j=0; j<i; j++){
        r = 0.0;
        r += (Molecule.centers[i].x-Molecule.centers[j].x)*
             (Molecule.centers[i].x-Molecule.centers[j].x);
        r += (Molecule.centers[i].y-Molecule.centers[j].y)*
             (Molecule.centers[i].y-Molecule.centers[j].y);
        r += (Molecule.centers[i].z-Molecule.centers[j].z)*
             (Molecule.centers[i].z-Molecule.centers[j].z);
        r = 1.0/sqrt(r);
        E += (Molecule.centers[i].Z_nuc*Molecule.centers[j].Z_nuc)*r;
      }
  Molecule.Enuc = E;

  return;
}
