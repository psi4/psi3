#ifndef IWL_H
#define IWL_H

#include <stdio.h>
#include <psio.h>

/*
** IWL.H
** Header file for Integrals With Labels Library
**
** David Sherrill
** Center for Computational Quantum Chemistry, UGA
**
*/

typedef short int Label;
typedef double Value;

struct iwlbuf {
   int itap;                   /* tape number for input file */
   psio_address bufpos;        /* current page/offset */
   int ints_per_buf;           /* integrals per buffer */
   int bufszc;                 /* buffer size in characters (bytes) */
   double cutoff;              /* cutoff value for writing */
   int lastbuf;                /* is this the last IWL buffer? 1=yes,0=no */
   int inbuf;                  /* how many ints in current buffer? */
   int idx;                    /* index of integral in current buffer */
   Label *labels;              /* pointer to where integral values begin */
   Value *values;              /* integral values */
   };

#define IWL_KEY_BUF "IWL Buffers"
#define IWL_KEY_EFZC "IWL Frozen Core Energy"
#define IWL_KEY_ONEL "IWL One-electron matrix elements"

#define IWL_INTS_PER_BUF 2980

extern void iwl_buf_fetch(struct iwlbuf *Buf);
extern int iwl_rdone(int itap, double *ints, double *e_fzc, int *ioff,
              int norbs, int nfzc, int nfzv, int erase,
              int printflg, FILE *outfile);
extern int iwl_rdone_all(int itap, int nbstri, double *onel_ints,
              double *efzc, int erase);
/*
extern void iwl_rdone(int itap, int nbstri, double *onel_ints, double *efzc,
      int erase);
*/
extern void iwl_wrtone(int itap, int ntri, double *onel_ints, double e_fzc);
extern void iwl_rdtwo(int itap, double *ints, int *ioff, int norbs,
                      int nfzc, int nfzv, int printflg, FILE *outfile);
extern void iwl_wrttwo(int itap, int nbfso, double *ints, int *ioff, 
      double toler, int printflg, FILE *outfile);
extern void sortbuf(struct iwlbuf *inbuf, struct iwlbuf *outbuf,
      double *ints, int fpq, int lpq, int *ioff, int *ioff2, 
      int nbfso, int elbert, int intermediate, int no_pq_perm, 
      int qdim, int add, int printflg, FILE *outfile); 
extern void iwl_buf_init(struct iwlbuf *Buf, int intape, double cutoff,
      int oldfile, int readflag);
extern int iwl_buf_rd(struct iwlbuf *Buf, int target_pq, double *ints,
      int *ioff_lt, int *ioff_rt, int mp2, int printflg, FILE *outfile);
extern int iwl_buf_rd_all(struct iwlbuf *Buf, double *ints,
      int *ioff_lt, int *ioff_rt, int mp2, int *ioff,
      int printflg, FILE *outfile);
extern int iwl_buf_rd_all_act(struct iwlbuf *Buf, double *ints,
      int *ioff_lt, int *ioff_rt, int mp2, int *ioff,
      int fstact, int lstact, int printflg, FILE *outfile);
extern int iwl_buf_rd_all_mp2r12a(struct iwlbuf *Buf, double *ints,
      int *ioff_lt, int *ioff_rt, int bra_ket_symm, int *ioff,
      int printflg, FILE *outfile);
extern void iwl_buf_wrt_all(struct iwlbuf *Buf, int nbfso, double *ints, 
      int *ioff, int printflg, FILE *outfile);
extern void iwl_buf_wrt(struct iwlbuf *Buf, int p, int q, int pq, int pqsym,
      double *arr, int rmax, int *active, int *ioff, int *orbsym, int *firsti, 
      int *lasti, int sortby_rs, int printflag, FILE *outfile);
extern void iwl_buf_wrt_mp2(struct iwlbuf *Buf, int p, int q, int pq,
      int pqsym, double **arr, int rsym, int *firstr, int *lastr, 
      int *firsts, int *lasts, int *occ, int *vir, int *ioff, 
      int printflag, FILE *outfile);
extern void iwl_buf_wrt_mp2r12a(struct iwlbuf *Buf, int p, int q, int pq,
      int pqsym, double **arr, int rsym, int *firstr, int *lastr, 
      int *firsts, int *lasts, int *occ, int bra_ket_symm, int *ioff, 
      int printflag, FILE *outfile);
extern void iwl_buf_flush(struct iwlbuf *Buf, int lastbuf);
extern void iwl_buf_close(struct iwlbuf *Buf, int keep);
extern void iwl_buf_toend(struct iwlbuf *Buf);
extern void iwl_buf_wrt_arr(struct iwlbuf *Buf, double *arr, int *p, int *q,
      int *r, int *s, int size);
extern void iwl_buf_wrt_arr_SI(struct iwlbuf *Buf, double *arr, 
      short int *p, short int *q, short int *r, short int *s, int size);
extern void iwl_buf_wrt_arr_SI_nocut(struct iwlbuf *Buf, double *arr,
      short int *p, short int *q, short int *r, short int *s, int size);
extern int iwl_buf_rd_arr(struct iwlbuf *Buf, int target_pq, double *ints,
      int *rlist, int *slist, int *size, int *ioff,
      int printflg, FILE *outfile);
extern int iwl_buf_rd_arr2(struct iwlbuf *Buf, double *ints, int *plist,
      int *qlist, int *rlist, int *slist, int *size, int *ioff,
      int printflg, FILE *outfile);
extern void iwl_buf_wrt_arr2(struct iwlbuf *Buf, double *arr, int p, int q, 
      int *rlist, int *slist, int size, int printflag, FILE *outfile);
extern void iwl_buf_wrt_mat(struct iwlbuf *Buf, int ptr, int qtr,
      double **mat, int rfirst, int rlast, int sfirst, int slast,
      int *reorder, int reorder_offset, int printflag, int *ioff,
      FILE *outfile);
extern void iwl_buf_wrt_val(struct iwlbuf *Buf, int p, int q, int r, int s,
                     double value, int printflag, FILE *outfile, int dirac);
extern void iwl_buf_wrt_val_SI(struct iwlbuf *Buf, short int p, short int q,
                     short int r, short int s, double value, int printflag,
                     FILE *outfile, int dirac);

#endif /* end IWL_H */
