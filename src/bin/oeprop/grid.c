#define EXTERN
#include "includes.h"
#include "globals.h"

extern void compute_grid_oeprops();
extern void compute_grid_dens();
extern void compute_grid_mos();


void compute_grid() {

  switch (grid) {
  case 1:
      compute_grid_oeprops();
      break;
      
  case 2:
  case 3:
  case 4:
      compute_grid_dens();
      break;

  case 5:
      compute_grid_mos();
      break;
  }

  return;
}
