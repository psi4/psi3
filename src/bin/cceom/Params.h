/* Input parameters for ccenergy */
struct Params {
    long int memory;
    int aobasis;
    int cachelev;
    int cachetype;
    int ref;
};

struct Eom_params {
  int max_iter;
  int vectors_per_root;
  int *rpi;
  double eval_tol;
  double residual_tol;
  int prop_root;
  int print_singles;
  double complex_tol;
  double schmidt_add_residual_tol;

  int max_iter_SS;
  int vectors_per_root_SS;
  int excitation_range;
  double residual_tol_SS;
};

