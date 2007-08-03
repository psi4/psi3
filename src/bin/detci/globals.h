/*! \file 
    \ingroup (DETCI)
    \brief Enter brief description of file here 
*/
/* nice stuff to extern or not to extern properly */
#ifdef EXTERN
# undef EXTERN
# define EXTERN extern
#else
# define EXTERN
#endif

EXTERN int errcod;
EXTERN struct calcinfo CalcInfo;
EXTERN struct params Parameters;
EXTERN int *ioff;
EXTERN FILE *infile, *outfile;
EXTERN char *psi_file_prefix;
EXTERN struct ci_blks CIblks;
EXTERN struct olsen_graph *AlphaG;
EXTERN struct olsen_graph *BetaG;
EXTERN struct graph_set *AlphaGraph;
EXTERN struct graph_set *BetaGraph;
EXTERN struct H_zero_block H0block;
EXTERN int ***OV;
EXTERN int **s1_contrib, **s2_contrib, **s3_contrib;
EXTERN double *tmp_ras_array;
EXTERN struct detci_timings detci_time;
