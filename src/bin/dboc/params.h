
#ifndef _psi3_bin_dboc_params_h_
#define _psi3_bin_dboc_params_h_

namespace PrintLevels {
  static const int print_intro = 1;
  static const int print_params = 1;
  static const int print_contrib = 2;
};

typedef struct {

  enum RefType { rhf=1, rohf=2, uhf=3};

  typedef struct {
    int index;
    double coeff;
  } Coord_t;

  char *label;
  char *wfn;
  RefType reftype;
  double delta;
  int disp_per_coord;
  int ncoord;
  Coord_t* coords;
  int nisotope;
  char** isotopes;

  int print_lvl;

} Params_t;

#endif
