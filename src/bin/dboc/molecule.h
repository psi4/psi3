
#ifndef _psi3_bin_dboc_molecule_h_
#define _psi3_bin_dboc_molecule_h_

typedef struct {
  int natom;
  double **geom;
  double *zvals;
  double *masses;
} Molecule_t;

#endif

