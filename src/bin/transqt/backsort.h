#ifndef BACKSORT_H
#define BACKSORT_H

void backsort_prep(int uhf);
void backsort(int first_tmp_file, double tolerance, int uhf);
void backsort_write(int i, int j, double **A, int kfirst, int klast,
		    int lfirst, int llast, int printflag, FILE *outfile,
		       struct iwlbuf *twopdm_out, int uhf);

#endif /* BACKSORT_H */
