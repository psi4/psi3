#ifndef DPD_H
#define DPD_H

#include <stdio.h>
#include <psio.h>
#include <libciomr.h>
#include <iwl.h>

struct dpdparams {
    int nirreps;      /* No. of irreps */
    int pqnum;
    int rsnum;
    int *rowtot;      /* Row dimension for each submatrix */
    int *coltot;      /* Column dimension for each submatrix */
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
};

struct dpdfile {
    char label[PSIO_KEYLEN];
    int filenum;
    psio_address *lfiles;
    struct dpdparams *params;
    double ***matrix;
};

struct dpdshift {
    int shift_type;
    int **rowtot;
    int **coltot;
    double ****matrix;
};

struct dpdbuf {
    int anti;
    struct dpdparams *params;
    struct dpdshift shift;
    struct dpdfile file;
    double ***matrix;
};

struct dpdtrans {
    double ***matrix;
    struct dpdshift shift;
    struct dpdbuf buf;
};

struct oe_dpdparams {
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
};

struct oe_dpdfile {
    char label[PSIO_KEYLEN];
    int filenum;
    psio_address *lfiles;
    struct oe_dpdparams *params;
    double ***matrix;
};

/* Useful for the generalized sorting function */
enum indices {pqrs, pqsr, prqs, prsq, psqr, psrq,
	      qprs, qpsr, qrps, qrsp, qspr, qsrp,
	      rqps, rqsp, rpqs, rpsq, rsqp, rspq,
	      sqrp, sqpr, srqp, srpq, spqr, sprq};

int dpd_init(int nirreps, int memory, int num_subspaces, ...);

int dpd_close(void);

int dpd_axpy(struct dpdbuf *BufX, struct dpdbuf *BufY, double alpha,
	     int print_flag, FILE *outfile);

int dpd_buf_close(struct dpdbuf *Buf);

int dpd_buf_dump(struct dpdbuf *DPDBuf, struct iwlbuf *IWLBuf,
		 int *prel, int *qrel, int *rrel, int *srel, 
		 int bk_pack, int swap23, int print_flag, FILE *outfile);

int dpd_buf_init(struct dpdbuf *Buf, int inputfile, int pqnum, int rsnum,
		 int file_pqnum, int file_rsnum, int anti, char *label,
		 int print_flag, FILE *outfile);

int dpd_buf_mat_irrep_close(struct dpdbuf *Buf, int irrep);

int dpd_buf_mat_irrep_close_block(struct dpdbuf *Buf, int irrep, int num_pq);

int dpd_buf_mat_irrep_init(struct dpdbuf *Buf, int irrep);

int dpd_buf_mat_irrep_init_block(struct dpdbuf *Buf, int irrep, int num_pq);

int dpd_buf_mat_irrep_print(struct dpdbuf *Buf, int irrep, FILE *outfile);

int dpd_buf_mat_irrep_rd(struct dpdbuf *Buf, int irrep, int print_flag,
			 FILE *outfile);

int dpd_buf_mat_irrep_rd_block(struct dpdbuf *Buf, int irrep, int start_pq,
			       int num_pq, int print_flag, FILE *outfile);

int dpd_buf_mat_irrep_row_close(struct dpdbuf *Buf, int irrep);

int dpd_buf_mat_irrep_row_init(struct dpdbuf *Buf, int irrep);

int dpd_buf_mat_irrep_row_rd(struct dpdbuf *Buf, int irrep, int pq,
			     int print_flag, FILE *outfile);

int dpd_buf_mat_irrep_row_wrt(struct dpdbuf *Buf, int irrep, int pq,
			      int print_flag, FILE *outfile);


int dpd_buf_mat_irrep_row_zero(struct dpdbuf *Buf, int irrep, int row);

int dpd_buf_mat_irrep_shift13(struct dpdbuf *Buf, int irrep, int print_flag,
                              FILE *outfile);

int dpd_buf_mat_irrep_shift31(struct dpdbuf *Buf, int irrep, int print_flag,
                              FILE *outfile);

int dpd_buf_mat_irrep_wrt(struct dpdbuf *Buf, int irrep, int print_flag,
			  FILE *outfile);

int dpd_buf_mat_irrep_wrt_block(struct dpdbuf *Buf, int irrep, int start_pq,
				int num_pq, int print_flag, FILE *outfile);

int dpd_buf_print(struct dpdbuf *Buf, FILE *outfile);

int dpd_buf_sort(struct dpdbuf *InBuf, int outfilenum, enum indices index,
		 int pqnum, int rsnum, char *label, int print_flag,
		 FILE *outfile);

int dpd_buf_symm(struct dpdbuf *Buf);

int dpd_buf_symm2(struct dpdbuf *Buf1, struct dpdbuf *Buf2);

double dpd_chksum(struct dpdbuf *BufX);

int dpd_contract111(struct oe_dpdfile *X, struct oe_dpdfile *Y, struct
		    oe_dpdfile *Z, int target_X, int target_Y, double
		    alpha, double beta, int print_flag, FILE *outfile);

int dpd_contract121(struct dpdbuf *X, struct oe_dpdfile *Y, struct
		    oe_dpdfile *Z, int trans_Y, int trans_Z, double
		    alpha, double beta, int print_flag, FILE *outfile);

int dpd_contract122(struct dpdbuf *X, struct dpdbuf *Y, struct oe_dpdfile *Z,
		    int target_X, int target_Y, double alpha, double beta,
		    int print_flag, FILE *outfile);

int dpd_contract212(struct oe_dpdfile *X, struct dpdbuf *Y, struct
		    dpdbuf *Z, int sum_X, int sum_Y, int trans_Z, double
		    alpha, double beta, int print_flag, FILE *outfile);

int dpd_contract221(struct dpdbuf *X, struct oe_dpdfile *Y, struct
		    dpdbuf *Z, int sum_X, int sum_Y, int trans_Z, double
		    alpha, double beta, int print_flag, FILE *outfile);

int dpd_contract222(struct dpdbuf *X, struct dpdbuf *Y, struct dpdbuf *Z,
		    int target_X, int target_Y, double alpha,
		    double beta, int print_flag, FILE *outfile);

int dpd_copy(struct dpdbuf *InBuf, int outfilenum, 
	     char *label, int print_flag, FILE *outfile);

int dpd_dirprd(struct dpdbuf *BufA, struct dpdbuf *BufB, 
	       int print_flag, FILE *outfile);

double dpd_dot(struct dpdbuf *BufA, struct dpdbuf *BufB, 
	       int print_flag, FILE *outfile);

int dpd_dot13(struct oe_dpdfile *T, struct dpdbuf *I, struct oe_dpdfile *Z,
	      int transt, int transz, double alpha, double beta,
	      int print_flag, FILE *outfile);

int dpd_dot14(struct oe_dpdfile *T, struct dpdbuf *I, struct oe_dpdfile *Z,
	      int transt, int transz, double alpha, double beta,
	      int print_flag, FILE *outfile);

int dpd_dot23(struct oe_dpdfile *T, struct dpdbuf *I, struct oe_dpdfile *Z,
	      int transt, int transz, double alpha, double beta,
	      int print_flag, FILE *outfile);

int dpd_dot24(struct oe_dpdfile *T, struct dpdbuf *I, struct oe_dpdfile *Z,
	      int transt, int transz, double alpha, double beta,
	      int print_flag, FILE *outfile);

void dpd_error(char *caller, FILE *outfile);

int dpd_file_close(struct dpdfile *File);

int dpd_file_init(struct dpdfile *File, int filenum, int pqnum,
                  int rsnum,  char *label, int print_flag,
                  FILE *outfile);

int dpd_file_mat_close(struct dpdfile *File);

int dpd_file_mat_init(struct dpdfile *File);

int dpd_file_mat_irrep_close(struct dpdfile *File, int irrep);

int dpd_file_mat_irrep_close_block(struct dpdfile *File, int irrep, int num_pq);

int dpd_file_mat_irrep_init(struct dpdfile *File, int irrep);

int dpd_file_mat_irrep_init_block(struct dpdfile *File, int irrep, int num_pq);

int dpd_file_mat_irrep_print(struct dpdfile *File, int irrep, FILE *outfile);

int dpd_file_mat_irrep_rd(struct dpdfile *File, int irrep, int print_flag,
			  FILE *outfile);

int dpd_file_mat_irrep_rd_block(struct dpdfile *File, int irrep, int start_pq,
				int num_pq, int print_flag, FILE *outfile);

int dpd_file_mat_irrep_row_close(struct dpdfile *File, int irrep);

int dpd_file_mat_irrep_row_init(struct dpdfile *File, int irrep);

int dpd_file_mat_irrep_row_rd(struct dpdfile *File, int irrep, int row);

int dpd_file_mat_irrep_row_wrt(struct dpdfile *File, int irrep, int row);

int dpd_file_mat_irrep_row_zero(struct dpdfile *File, int irrep, int row);

int dpd_file_mat_irrep_wrt(struct dpdfile *File, int irrep, int print_flag,
                           FILE *outfile);

int dpd_file_mat_irrep_wrt_block(struct dpdfile *File, int irrep, int start_pq,
				 int num_pq, int print_flag, FILE *outfile);

int dpd_file_mat_wrt(struct dpdfile *File, int print_flag, FILE *outfile);

int dpd_file_print(struct dpdfile *File, FILE *outfile);

int dpd_mat_irrep_print(double **matrix, struct dpdparams *Params,
                        int irrep, FILE *outfile);

int dpd_oe_axpy(struct oe_dpdfile *FileA, struct oe_dpdfile *FileB, 
		double alpha, int transa, int print_flag, FILE *outfile);

double dpd_oe_chksum(struct oe_dpdfile *BufX);

int dpd_oe_copy(struct oe_dpdfile *InFile, int outfilenum, 
		char *label, int print_flag, FILE *outfile);

int dpd_oe_dirprd(struct oe_dpdfile *FileA, struct oe_dpdfile *FileB, 
		  int print_flag, FILE *outfile);

double dpd_oe_dot(struct oe_dpdfile *FileA, struct oe_dpdfile *FileB, 
		  int print_flag, FILE *outfile);

int dpd_oe_file_close(struct oe_dpdfile *File);

int dpd_oe_file_init(struct oe_dpdfile *File, int filenum, int pnum, int qnum,
		     char *label, int print_flag, FILE *outfile);

int dpd_oe_file_mat_close(struct oe_dpdfile *File);

int dpd_oe_file_mat_init(struct oe_dpdfile *File);

int dpd_oe_file_mat_irrep_close(struct oe_dpdfile *File, int irrep);

int dpd_oe_file_mat_irrep_init(struct oe_dpdfile *File, int irrep);

int dpd_oe_file_mat_irrep_print(struct oe_dpdfile *File, int irrep, FILE *outfile);

int dpd_oe_file_mat_irrep_rd(struct oe_dpdfile *File, int irrep, int print_flag,
			     FILE *outfile);

int dpd_oe_file_mat_irrep_wrt(struct oe_dpdfile *File, int irrep, int print_flag,
			      FILE *outfile);

int dpd_oe_file_mat_rd(struct oe_dpdfile *File, int print_flag, FILE *outfile);

int dpd_oe_file_mat_wrt(struct oe_dpdfile *File, int print_flag, FILE *outfile);

int dpd_oe_file_print(struct oe_dpdfile *File, FILE *outfile);

int dpd_oe_mat_irrep_print(double **matrix, struct oe_dpdparams *Params,
			   int irrep, FILE *outfile);

int dpd_oe_params_print(struct oe_dpdparams *Params, FILE *outfile);

int dpd_oe_scm(struct oe_dpdfile *InFile, double alpha, int print_flag,
	       FILE *outfile);

double dpd_oe_trace(struct oe_dpdfile *InFile, int print_flag, FILE *outfile);

int dpd_params_print(struct dpdparams *Params, FILE *outfile);

int dpd_scm(struct dpdbuf *InBuf, double alpha, int print_flag, FILE *outfile);

int dpd_swapbk(struct dpdbuf *InBuf, int outfilenum, int pqnum, int rsnum,
	       char *label, int print_flag, FILE *outfile);

int dpd_swap12(struct dpdbuf *InBuf, int outfilenum, int pqnum, int rsnum,
	       char *label, int print_flag, FILE *outfile);

int dpd_swap13(struct dpdbuf *InBuf, int outfilenum, int pqnum, int rsnum,
	       char *label, int print_flag, FILE *outfile);

int dpd_swap14(struct dpdbuf *InBuf, int outfilenum, int pqnum, int rsnum,
	       char *label, int print_flag, FILE *outfile);

int dpd_swap23(struct dpdbuf *InBuf, int outfilenum, int pqnum, int rsnum,
	       char *label, int print_flag, FILE *outfile);

int dpd_swap24(struct dpdbuf *InBuf, int outfilenum, int pqnum, int rsnum,
	       char *label, int print_flag, FILE *outfile);

int dpd_swap34(struct dpdbuf *InBuf, int outfilenum, int pqnum, int rsnum,
	       char *label, int print_flag, FILE *outfile);

int dpd_trans_close(struct dpdtrans *Trans);

int dpd_trans_init(struct dpdtrans *Trans, struct dpdbuf *Buf);

int dpd_trans_mat_irrep_close(struct dpdtrans *Trans, int irrep);

int dpd_trans_mat_irrep_init(struct dpdtrans *Trans, int irrep);

int dpd_trans_mat_irrep_rd(struct dpdtrans *Trans, int irrep);

int dpd_trans_mat_irrep_shift13(struct dpdtrans *Trans, int irrep,
				int print_flag, FILE *outfile);

int dpd_trans_mat_irrep_shift31(struct dpdtrans *Trans, int irrep,
				int print_flag, FILE *outfile);

int dpd_trans_mat_irrep_wrt(struct dpdtrans *Trans, int irrep);

#endif /* DPD_H */

