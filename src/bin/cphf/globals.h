
/* Global variables */
#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

EXTERN FILE *infile, *outfile;
EXTERN int *ioff;
#define IOFF_MAX 32641
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

/* setup.c */
EXTERN int natom, nmo, nso, nao, nirreps, ndocc, nuocc;
EXTERN int ntri, num_ai, num_ij, ntei, noei, noei_ao;
EXTERN int *orbspi, *clsdpi, *uoccpi;
EXTERN double *evals, *zvals, *ints, **scf, **usotao, **geom;
EXTERN int *first, *last, *ofirst, *olast, *vfirst, *vlast;

EXTERN int print_lvl;
