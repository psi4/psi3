#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<libciomr.h>
#include<file30.h>

#include<libint.h>
#include"defines.h"
#define EXTERN
#include"global.h"

/*-------------------------------
  Explicit function declarations
 -------------------------------*/
static void get_geometry(void);
static void compute_enuc(void);

void init_molecule()
{
  Molecule.label = file30_rd_label();
  Molecule.num_atoms = file30_rd_natom();
/* Molecule.centers = */ get_geometry();

  return;
}


void cleanup_molecule()
{
  free(Molecule.centers);

  return;
}

void get_geometry()
{
   int i;
   double *Z;  /* nuclear charges */
   double **g; /* cartesian geometry */
   
   Molecule.centers = (struct coordinates *)malloc(sizeof(struct coordinates)*Molecule.num_atoms);
   g = file30_rd_geom();
   Z = file30_rd_zvals();

   /*--- move it into the appropriate struct form ---*/
   for (i=0; i<Molecule.num_atoms; i++){
      Molecule.centers[i].x = g[i][0];
      Molecule.centers[i].y = g[i][1];
      Molecule.centers[i].z = g[i][2];
      Molecule.centers[i].Z_nuc = Z[i];
   }
   free_block(g);
   free(Z);

   return;
}


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
