#ifdef EXTERN
# undef EXTERN
# define EXTERN extern
extern int MAX_LINELENGTH;
#else 
# define EXTERN
 int MAX_LINELENGTH = 133;
#endif



EXTERN FILE *infile, *outfile;

EXTERN int num_atoms;
EXTERN double **cart_geom;
EXTERN double **cart_grad;
EXTERN double *atomic_nums;

EXTERN void punt(char *mess);
EXTERN void read_file11();

/*this needs to be in C*/
extern "C" char *gprgid(); 
