#include <stdio.h>
#include <stdlib.h>
extern "C" {
#include <libciomr/libciomr.h>
}
#include "molecule.h"
#include "params.h"

extern Molecule_t Molecule;
extern Params_t Params;

static void append_geom(FILE *geomdat, double **geom, int disp);

void setup_geoms()
{
  const int disp_per_coord = Params.disp_per_coord;
  const double delta = Params.delta;
  int d;
  int disp, atom, xyz;
  FILE *geometry;

  /*--- Open geom.dat for writing ---*/
  ffile(&geometry, "geom.dat", 0);

  /*--- make a local copy of reference geometry ---*/
  double **geom_copy = block_matrix(Molecule.natom,3);
  for(atom=0; atom<Molecule.natom; atom++)
    for(xyz=0; xyz<3; xyz++)
      geom_copy[atom][xyz] = Molecule.geom[atom][xyz];

  
  for(atom=0,disp=1; atom<Molecule.natom; atom++)
    for(xyz=0; xyz<3; xyz++)
      for(d=-1; d<=1; d+=2,disp++) {
	geom_copy[atom][xyz] += d*delta;
	append_geom(geometry,geom_copy,disp);
	geom_copy[atom][xyz] -= d*delta;
      }
  
  free_block(geom_copy);
  fclose(geometry);

}


void append_geom(FILE *geomdat, double **geom, int disp)
{
  fprintf(geomdat,"%%% DBOC cartesian displacement %d\n",disp);
  fprintf(geomdat,"geometry%d = (\n",disp);
  for(int atom=0; atom<Molecule.natom; atom++)
    fprintf(geomdat,"  (%4.2lf %20.15lf %20.15lf %20.15lf)\n",
	    Molecule.zvals[atom],geom[atom][0],geom[atom][1],geom[atom][2]);
}
