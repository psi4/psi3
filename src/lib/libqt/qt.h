/*
** qt.h
**
** Header file for the Quantum Trio Library 
**
** David Sherrill 1994
**
** Modifications by Daniel Crawford 1996, 1997
*/

int mat_in(FILE *fp, double **array, int width, int max_length, int *stat);
void fill_sym_matrix(double **A, int size);
double combinations(int n, int k);
double factorial(int n);
void schmidt(double **A, int rows, int cols, FILE *outfile);
int schmidt_add(double **A, int rows, int cols, double *v);
void normalize(double **A, int rows, int cols);
double invert_matrix(double **a, double **y, int N, FILE *outfile);
void solve_2x2_pep(double **H, double S, double *evals, double **evecs);
void reorder_qt(int *docc_in, int *socc_in, int *frozen_docc_in,
      int *frozen_uocc_in, int *order, int *orbs_per_irrep, int nirreps);
void reorder_qt_uhf(int *docc, int *socc, int *frozen_docc, 
		    int *frozen_uocc, int *order_alpha, int *order_beta,
		    int *orbspi, int nirreps);
void reorder_ras(int *docc_in, int *socc_in, int *frozen_docc_in,
      int *frozen_uocc_in, int *order, int *orbs_per_irrep,
      int *ras1, int *ras2, int *ras3, int *ras4, int do_ras4, int nirreps);
void reorder_ras2(int *docc_in, int *socc_in, int *frozen_docc_in, 
      int *frozen_uocc_in, int *order, int *orbs_per_irrep, 
      int *ras1, int *ras2, int *ras3, int *ras4, int parsed_ras1,
      int parsed_ras2, int do_ras4, int nirreps);
int ras_set(int nirreps, int nbfso, int freeze_core, int *orbspi,
      int *docc, int *socc, int *frdocc, int *fruocc,
      int **ras_opi, int *order, int ras_type);
int ras_set2(int nirreps, int nbfso, int delete_fzdocc,
      int delete_restrdocc, int *orbspi,
      int *docc, int *socc, int *frdocc, int *fruocc,
      int *restrdocc, int *restruocc, int **ras_opi, int *order,
      int ras_type);
void newmm_rking(double **A, int transa, double **B, int transb, double **C,
           int num_rows, int num_links, int num_cols, 
           double alpha, double beta);
double dot_block(double **A, double **B, int rows, int cols, double alpha);
void dirprd_block(double **A, double **B, int rows, int cols);
int pople(double **A, double *x, int dimen, int num_vecs, double tolerance,
           FILE *outfile, int print_lvl);
void mat_print(double **A, int rows, int cols, FILE *outfile);
double eri(unsigned int l1, unsigned int m1, unsigned int n1, double alpha1, double A[3],
	   unsigned int l2, unsigned int m2, unsigned int n2, double alpha2, double B[3],
	   unsigned int l3, unsigned int m3, unsigned int n3, double alpha3, double C[3],
	   unsigned int l4, unsigned int m4, unsigned int n4, double alpha4, double D[3],
	   int norm_flag);
double norm_const(unsigned int l1, unsigned int m1, unsigned int n1, double alpha1, double A[3]);

extern void C_DAXPY(int length, double a, double *x, int inc_x,
             double *y, int inc_y);
extern void C_DCOPY(int length, double *x, int inc_x,
             double *y, int inc_y);
extern void C_DGEMM(char transa, char transb, int m, int n, int k,
                    double alpha, double *A, int nca, double *B, int ncb,
                    double beta, double *C, int ncc);
extern void C_DROT(int ntot, double *x, int incx, double *y, int incy,
                   double costheta, double sintheta);
extern void C_DSCAL(int len, double alpha, double *vec, int inc);
extern void C_DGEMV(char transa, int m, int n, double alpha, double *A, 
		    int nca, double *X, int inc_x, double beta, double *Y,
		    int inc_y);
extern double C_DDOT(int n, double *X, int inc_x, double *Y, int inc_y);
extern int C_DGETRF(int nrow, int ncol, double *a, int lda, int *ipiv);
void timer_init(void);
void timer_done(void);
void timer_on(char *key);
void timer_off(char *key);

void filter(double *input, double *output, int *ioff, int norbs, int nfzc, int nfzv);

int C_DGEEV(int n, double **a, int lda,
  double *wr, double *wi, double **vl, int ldvl, double **vr,
  int ldvr, double *work, int lwork, int info);

int C_DGESV(int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb);
int C_DGETRI(int n, double *a, int lda, int *ipiv, double *work, int lwork);
int C_DGESVD(char jobu, char jobvt, int m, int n, double *A, int lda, double *s, int lds,
	     double *u, int ldu, double *vt, int ldvt, double *work, int lwork);
int C_DSYEV(char jobz, char uplo, int n, double *A, int lda, double *w, double *work, int lwork);

void print_block(double *, int, int, FILE *);

void sort(double *A, double **B, int n);

int david(double **A, int N, int M, double *eps, double **v, double cutoff, 
	  int print);
void init_in_out(int nstrings, char *fnames[]);

