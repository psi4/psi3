
/* $Id$ */
/* $Log$
 * Revision 1.1  2000/02/04 22:50:46  evaleev
 * Initial revision
 *
/* Revision 1.1  1991/06/15 22:45:47  seidl
/* Initial revision
/* */


#define MAX_BASIS 200
#define MAX_STRING 512
#ifdef ULTRIX
# define MAX_WORDS 1875000
#elif defined AIX
# define MAX_WORDS 3600000
#else
# define MAX_WORDS 1000000
#endif

#ifdef EXTERN
# undef EXTERN
# define EXTERN extern
#else
# define EXTERN
#endif

EXTERN FILE *infile, *outfile;
EXTERN FILE *file11,*file15,*file17,*file18;

EXTERN int print;   
EXTERN int dipole;
EXTERN int conv;

EXTERN int n_so_typs;
EXTERN int nbasis;
EXTERN int nbfao;
EXTERN int nbfso;
EXTERN int nbatri;
EXTERN int nbstri;
EXTERN int ntri;
EXTERN int maxbuf;
EXTERN int nsect;
EXTERN int natom;
EXTERN int natom3;
EXTERN int natom3x;
EXTERN int iopen;
EXTERN int twocon;
EXTERN int ntypes;
EXTERN int ntype1;
EXTERN int nc;
EXTERN int no;
EXTERN int nocc;
EXTERN int nind;
EXTERN int ndep;

EXTERN int itap40,itap37;              
EXTERN int itap42,itap43,itap44;
EXTERN int work,work2;

EXTERN int ioff[1024];
EXTERN int sa_loc[MAX_BASIS];
EXTERN int ha_loc[MAX_BASIS];
EXTERN int fa_loc[MAX_BASIS];
EXTERN int ea_loc[MAX_BASIS];
EXTERN int ba_loc[MAX_BASIS];
EXTERN int ua_loc[MAX_BASIS];
EXTERN unsigned short motyp[MAX_BASIS];
EXTERN unsigned short nstart[10],nend[10],nsorb[10];
EXTERN int *block_locs;

EXTERN double enuc,escf;

EXTERN double *e_vals,*occ_num;
EXTERN double **e_vecs_so,**e_vecs_ao,**ao_to_so;
EXTERN double **alpa,**beta;
EXTERN double focc[10];
EXTERN double **so_vecs_k,**so_vecs_l;

/* tcscf stuff */
EXTERN double c1sq,c2sq,c12p,c1,c2,elec;
EXTERN double h11,h22,h12,h21;
EXTERN double dpm[3],dpn[3],h11f[3],h22f[3];
EXTERN double *ha11,*ha22,*ha12,*e1t;
EXTERN double *c1a,*c2a;
EXTERN double a22[2][2];
EXTERN struct off_diags {
          double one;
          double two;
          } *baf,*a12;

/* end tcscf   */

EXTERN struct ind_pairs {
          unsigned short ii;
          unsigned short jj;
          unsigned int ij;
          unsigned short it;
          unsigned short jt;
          } *indep;

EXTERN struct o_pkints {
          unsigned int ij;
          unsigned int kl;
          double pval;
          double kval;
          } *o_pkbuf;

EXTERN struct c_pkints {
          unsigned int ij;
          unsigned int kl;
          double pval;
          } *c_pkbuf;
