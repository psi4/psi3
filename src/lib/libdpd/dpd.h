#ifndef DPD_H
#define DPD_H

#include <stdio.h>
#include <psio.h>
#include <libciomr.h>
#include <iwl.h>

typedef struct {
    int nirreps;      /* No. of irreps */
    int pqnum;        /* Pair number for the row indices */
    int rsnum;        /* Pair number for the column indices */
    int *rowtot;      /* Row dimension for each irrep */
    int *coltot;      /* Column dimension for each irrep */
    int **rowidx;     /* Row index lookup array */
    int **colidx;     /* Column index lookup array */
    int ***roworb;    /* Row index -> orbital index lookup array */
    int ***colorb;    /* Column index -> orbital index lookup array */
    int *ppi;         /* Number of p indices per irrep */
    int *qpi;         /* Number of q indices per irrep */
    int *rpi;         /* Number of r indices per irrep */
    int *spi;         /* Number of s indices per irrep */
    int *poff;        /* Orbital offset for p */
    int *qoff;        /* Orbital offset for q */
    int *roff;        /* Orbital offset for r */
    int *soff;        /* Orbital offset for s */
    int *psym;        /* Orbital symmetry for index p */
    int *qsym;        /* Orbital symmetry for index q */
    int *rsym;        /* Orbital symmetry for index r */
    int *ssym;        /* Orbital symmetry for index s */
    int perm_pq;      /* Can p and q be permuted? */
    int perm_rs;      /* Can r and s be permuted? */
    int peq;          /* Can p and q be equal? */
    int res;          /* Can r and s be equal? */
} dpdparams4;

typedef struct {
    char label[PSIO_KEYLEN];
    int filenum;
    int my_irrep;     /* Total irrep of this quantity */
    psio_address *lfiles;   /* File address for each submatrix by ROW irrep */
    dpdparams4 *params;
    int incore;
    double ***matrix;
} dpdfile4;

typedef struct {
    int shift_type;
    int **rowtot;
    int **coltot;
    double ****matrix;
} dpdshift4;

typedef struct {
    int anti;         /* Is this buffer antisymmetric? */
    dpdparams4 *params;
    dpdfile4 file;
    dpdshift4 shift;
    double ***matrix;
} dpdbuf4;

typedef struct {
    double ***matrix;
    dpdshift4 shift;
    dpdbuf4 buf;
} dpdtrans4;

typedef struct {
    int nirreps;      /* No. of irreps */
    int pnum;
    int qnum;
    int *rowtot;      /* Row dimension for each submatrix */
    int *coltot;      /* Column dimension for each submatrix */
    int *rowidx;      /* Row index lookup array */
    int *colidx;      /* Column index lookup array */
    int **roworb;     /* Row index -> orbital index lookup array */
    int **colorb;     /* Column index -> orbital index lookup array */
    int *ppi;         /* Number of p indices per irrep */
    int *qpi;         /* Number of q indices per irrep */
    int *poff;        /* Orbital offset for p */
    int *qoff;        /* Orbital offset for q */
    int *psym;        /* Orbital symmetry for index p */
    int *qsym;        /* Orbital symmetry for index q */
} dpdparams2;

typedef struct {
    char label[PSIO_KEYLEN];
    int filenum;
    int my_irrep;
    psio_address *lfiles;
    dpdparams2 *params;
    int incore;
    double ***matrix;
} dpdfile2;

/* DPD File4 Cache entries */
struct dpd_file4_cache_entry {
    int filenum;
    int irrep;
    int pqnum;
    int rsnum;
    char label[PSIO_KEYLEN];
    double ***matrix;
    int size;
    struct dpd_file4_cache_entry *next;
    struct dpd_file4_cache_entry *last;
};

/* DPD File2 Cache entries */
struct dpd_file2_cache_entry {
    int filenum;
    int irrep;
    int pnum;
    int qnum;
    char label[PSIO_KEYLEN];
    double ***matrix;
    int size;
    struct dpd_file2_cache_entry *next;
    struct dpd_file2_cache_entry *last;
};

/* DPD global parameter set */
typedef struct {
    int nirreps;
    int memory;
    int num_subspaces;
    int num_pairs;
    int *numorbs;
    int **orboff;
    int **pairtot;
    int **orbspi;
    int **orbsym;
    int **orbidx2;
    int ***pairidx;
    int ***orbs2;
    int ****pairorb;
    dpdparams2 **params2;
    dpdparams4 **params4;
    struct dpd_file2_cache_entry *file2_cache;
    struct dpd_file4_cache_entry *file4_cache;
    int *cachefiles;
    int **cachelist;
} dpd_data;

/* Useful for the generalized 4-index sorting function */
enum indices {pqrs, pqsr, prqs, prsq, psqr, psrq,
	      qprs, qpsr, qrps, qrsp, qspr, qsrp,
	      rqps, rqsp, rpqs, rpsq, rsqp, rspq,
	      sqrp, sqpr, srqp, srpq, spqr, sprq};

int dpd_init(int dpd_num, int nirreps, int memory, int *cachefiles,
             int **cachelist, int num_subspaces, ...);
int dpd_close(int dpd_num);
int dpd_set_default(int dpd_num);

void dpd_error(char *caller, FILE *outfile);

int dpd_contract222(dpdfile2 *X, dpdfile2 *Y, dpdfile2 *Z, int target_X,
		    int target_Y, double alpha, double beta);
int dpd_contract442(dpdbuf4 *X, dpdbuf4 *Y, dpdfile2 *Z, int target_X,
		    int target_Y, double alpha, double beta);
int dpd_contract422(dpdbuf4 *X, dpdfile2 *Y, dpdfile2 *Z, int trans_Y,
		    int trans_Z, double alpha, double beta);
int dpd_contract244(dpdfile2 *X, dpdbuf4 *Y, dpdbuf4 *Z, int sum_X, int sum_Y,
		    int trans_Z, double alpha, double beta);
int dpd_contract424(dpdbuf4 *X, dpdfile2 *Y, dpdbuf4 *Z, int sum_X,
		    int sum_Y, int trans_Z, double alpha, double beta);
int dpd_contract444(dpdbuf4 *X, dpdbuf4 *Y, dpdbuf4 *Z,
                    int target_X, int target_Y, double alpha, double beta);

/* Need to consolidate these routines into one general function */
int dpd_dot23(dpdfile2 *T, dpdbuf4 *I, dpdfile2 *Z,
	      int transt, int transz, double alpha, double beta);
int dpd_dot24(dpdfile2 *T, dpdbuf4 *I, dpdfile2 *Z,
	      int transt, int transz, double alpha, double beta);
int dpd_dot13(dpdfile2 *T, dpdbuf4 *I, dpdfile2 *Z,
	      int transt, int transz, double alpha, double beta);
int dpd_dot14(dpdfile2 *T, dpdbuf4 *I, dpdfile2 *Z,
	      int transt, int transz, double alpha, double beta);

int dpd_file2_init(dpdfile2 *File, int filenum, int irrep, int pnum,
		   int qnum, char *label);
int dpd_file2_close(dpdfile2 *File);
int dpd_file2_mat_init(dpdfile2 *File);
int dpd_file2_mat_close(dpdfile2 *File);
int dpd_file2_mat_rd(dpdfile2 *File);
int dpd_file2_mat_wrt(dpdfile2 *File);
int dpd_file2_print(dpdfile2 *File, FILE *outfile);
int dpd_file2_copy(dpdfile2 *InFile, int outfilenum, char *label);
int dpd_file2_dirprd(dpdfile2 *FileA, dpdfile2 *FileB);
double dpd_file2_dot(dpdfile2 *FileA, dpdfile2 *FileB);
int dpd_file2_scm(dpdfile2 *InFile, double alpha);
double dpd_file2_norm(dpdfile2 *BufX);
double dpd_file2_trace(dpdfile2 *InFile);
int dpd_file2_axpy(dpdfile2 *FileA, dpdfile2 *FileB,
                   double alpha, int transA);

int dpd_file4_init(dpdfile4 *File, int filenum, int irrep, int pqnum,
		   int rsnum,  char *label);
int dpd_file4_init_nocache(dpdfile4 *File, int filenum, int irrep, int pqnum,
		   int rsnum,  char *label);
int dpd_file4_close(dpdfile4 *File);
int dpd_file4_mat_irrep_init(dpdfile4 *File, int irrep);
int dpd_file4_mat_irrep_close(dpdfile4 *File, int irrep);
int dpd_file4_mat_irrep_rd(dpdfile4 *File, int irrep);
int dpd_file4_mat_irrep_wrt(dpdfile4 *File, int irrep);
int dpd_file4_mat_irrep_row_init(dpdfile4 *File, int irrep);
int dpd_file4_mat_irrep_row_close(dpdfile4 *File, int irrep);
int dpd_file4_mat_irrep_row_rd(dpdfile4 *File, int irrep, int row);
int dpd_file4_mat_irrep_row_wrt(dpdfile4 *File, int irrep, int row);
int dpd_file4_mat_irrep_row_zero(dpdfile4 *File, int irrep, int row);
int dpd_file4_print(dpdfile4 *File, FILE *outfile);
int dpd_file4_mat_irrep_rd_block(dpdfile4 *File, int irrep, int start_pq,
				int num_pq);
int dpd_file4_mat_irrep_wrt_block(dpdfile4 *File, int irrep, int start_pq,
				 int num_pq);

int dpd_buf4_init(dpdbuf4 *Buf, int inputfile, int irrep, int pqnum, int rsnum,
		 int file_pqnum, int file_rsnum, int anti, char *label);
int dpd_buf4_close(dpdbuf4 *Buf);
int dpd_buf4_mat_irrep_init(dpdbuf4 *Buf, int irrep);
int dpd_buf4_mat_irrep_close(dpdbuf4 *Buf, int irrep);
int dpd_buf4_mat_irrep_rd(dpdbuf4 *Buf, int irrep);
int dpd_buf4_print(dpdbuf4 *Buf, FILE *outfile);
int dpd_buf4_copy(dpdbuf4 *InBuf, int outfilenum, char *label);
int dpd_buf4_sort(dpdbuf4 *InBuf, int outfilenum, enum indices index,
		  int pqnum, int rsnum, char *label);
int dpd_buf4_sort_ooc(dpdbuf4 *InBuf, int outfilenum, enum indices index,
		      int pqnum, int rsnum, char *label);
int dpd_buf4_axpy(dpdbuf4 *BufX, dpdbuf4 *BufY, double alpha);
int dpd_buf4_dirprd(dpdbuf4 *BufA, dpdbuf4 *BufB);
double dpd_buf4_dot(dpdbuf4 *BufA, dpdbuf4 *BufB);
double dpd_buf4_norm(dpdbuf4 *BufX);
int dpd_buf4_scm(dpdbuf4 *InBuf, double alpha);
int dpd_buf4_scmcopy(dpdbuf4 *InBuf, int outfilenum, char *label, 
                     double alpha);
int dpd_buf4_symm(dpdbuf4 *Buf);
int dpd_buf4_symm2(dpdbuf4 *Buf1, dpdbuf4 *Buf2);
int dpd_buf4_mat_irrep_shift13(dpdbuf4 *Buf, int irrep);
int dpd_buf4_mat_irrep_shift31(dpdbuf4 *Buf, int irrep);
int dpd_buf4_mat_irrep_row_init(dpdbuf4 *Buf, int irrep);
int dpd_buf4_mat_irrep_row_close(dpdbuf4 *Buf, int irrep);
int dpd_buf4_mat_irrep_row_zero(dpdbuf4 *Buf, int irrep, int row);
int dpd_buf4_mat_irrep_row_rd(dpdbuf4 *Buf, int irrep, int pq);
int dpd_buf4_mat_irrep_row_wrt(dpdbuf4 *Buf, int irrep, int pq);
int dpd_buf4_mat_irrep_init_block(dpdbuf4 *Buf, int irrep, int num_pq);
int dpd_buf4_mat_irrep_close_block(dpdbuf4 *Buf, int irrep, int num_pq);
int dpd_buf4_mat_irrep_rd_block(dpdbuf4 *Buf, int irrep, int start_pq,
				int num_pq);
int dpd_buf4_mat_irrep_wrt_block(dpdbuf4 *Buf, int irrep, int start_pq,
				int num_pq);
int dpd_buf4_dump(dpdbuf4 *DPDBuf, struct iwlbuf *IWLBuf,
                  int *prel, int *qrel, int *rrel, int *srel,
                  int bk_pack, int swap23);

int dpd_trans4_init(dpdtrans4 *Trans, dpdbuf4 *Buf);
int dpd_trans4_close(dpdtrans4 *Trans);
int dpd_trans4_mat_irrep_init(dpdtrans4 *Trans, int irrep);
int dpd_trans4_mat_irrep_close(dpdtrans4 *Trans, int irrep);
int dpd_trans4_mat_irrep_rd(dpdtrans4 *Trans, int irrep);
int dpd_trans4_mat_irrep_wrt(dpdtrans4 *Trans, int irrep);
int dpd_trans4_mat_irrep_shift13(dpdtrans4 *Trans, int irrep);
int dpd_trans4_mat_irrep_shift31(dpdtrans4 *Trans, int irrep);

int dpd_4mat_irrep_print(double **matrix, dpdparams4 *Params,
			 int irrep, int my_irrep, FILE *outfile);

void dpd_file2_cache_init(void);
void dpd_file2_cache_close(void);
void dpd_file2_cache_print(FILE *outfile);
struct dpd_file2_cache_entry
 *dpd_file2_cache_scan(int filenum, int irrep, int pnum, int qnum, char *label);
struct dpd_file2_cache_entry *dpd_file2_cache_last(void);
int dpd_file2_cache_add(dpdfile2 *File);
int dpd_file2_cache_del(dpdfile2 *File);

void dpd_file4_cache_init(void);
void dpd_file4_cache_close(void);
void dpd_file4_cache_print(FILE *outfile);
struct dpd_file4_cache_entry
 *dpd_file4_cache_scan(int filenum, int irrep, int pqnum, int rsnum, char *label);
struct dpd_file4_cache_entry *dpd_file4_cache_last(void);
int dpd_file4_cache_add(dpdfile4 *File);
int dpd_file4_cache_del(dpdfile4 *File);

#endif /* DPD_H */

