#ifdef EXTERN
# undef EXTERN
# define EXTERN extern
#else
# define EXTERN
#endif

#define OPT_AUX (13)
#define SQR(A) ((A)*(A))
/*--- These are just internal type specifiers ---*/
#define STRE_TYPE (0)
#define BEND_TYPE (1)
#define TORS_TYPE (2)
#define OUT_TYPE (3)
#define NUM_INTCO_TYPES (4)
#define PRINT_TO_GEOM (113)
#define PRINT_TO_30 (114)
/* Limits to hard-wired arrays */
#define MAX_SALCS (1000)
#define MAX_ATOMS (300)
#define MAX_LINELENGTH (133)
#define MAX_SALC_LENGTH (1000)
#define MAX(I,J) ((I>J) ? I : J)
#define MIN(I,J) ((I>J) ? J : I)
#define EVAL_TOL (1.0E-14)                                     /* tolerance for eigenvalues (used in sq_rsp() and irrep() ) */
#define REDUNDANT_EVAL_TOL (1.0E-7)
#define SPANNED_IRREP_TOL (0.05)                               /* if character greater than this, irrep projected and kept */
#define LABEL_LENGTH (4) // for point group and irrep labels
/* step size limits */
#define STEP_LIMIT (0.1)                                       /* max step size if di coord has small value */
#define STEP_PERCENT (0.1)                                     /* if di coord large valued, max percentage allowed for step */

EXTERN FILE *fp_input, *fp_intco, *outfile, *fp_fconst, *fp_opt_aux;
EXTERN void cross_product(double *u,double *v,double *out);
EXTERN void scalar_mult(double a, double *vect, int dim);
EXTERN void scalar_div(double a, double *vect);
EXTERN double **symm_matrix_invert(double **A, int rows, int cols, int redundant);
EXTERN int div_int(int big, int little);
EXTERN void print_mat2(double **matrix, int rows, int cols, FILE *of);
EXTERN int *ops_in_class;
EXTERN int num_irreps, *irr;
EXTERN int num_nonzero;      /* number of non-redundant di coordinates (eigenvectors of G with nonzero eigenvalues) */
EXTERN char ptgrp[4];        /*molecular point group*/
EXTERN int *number_internals;
EXTERN void open_PSIF(void);
EXTERN void close_PSIF(void);
EXTERN double energy_chkpt(void);

/* C functions we may want to use (or need to keep as C functions so they can be called from C code) */
extern "C" void punt(char *message);
extern "C" char *gprgid();

struct OPTInfo {
  double *masses;

/* print options */
  int print_simples;
  int print_params;
  int print_delocalize;
  int print_symmetry;

/* optimization parameters */
  int optimize;
  int redundant;
  int delocalize;
  int bfgs;
  int dertype;
  int numerical_dertype;
  int iteration;
  int micro_iteration;
  double conv; /* MAX force */
  double ev_tol;
  double scale_connectivity;
  double edisp;
  int mix_types;

/* Back-transformation parameters */
  int bt_max_iter;
  double bt_dq_conv;
  double bt_dx_conv;

/* Obscure limits in intco evaluation */
  double cos_tors_near_1_tol;
  double cos_tors_near_neg1_tol;
  double sin_phi_denominator_tol;
};

struct SYMInfo {
  char *symmetry;
  int num_irreps;
  int **ct;
  char **irrep_lbls;
  char **clean_irrep_lbls;
  char **op_lbls;
  int **ict;
  int **ict_ops;
  int **ict_ops_sign;
};

EXTERN struct OPTInfo optinfo;
EXTERN struct SYMInfo syminfo;









