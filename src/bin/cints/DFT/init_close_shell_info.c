#include<stdio.h>
#include<stdlib.h>
#include <ip_libv1.h>
#include <libciomr.h>

#include"defines.h"
#define EXTERN
#include"global.h"

struct close_shell_info_s init_close_shell_info(void){
    
    struct close_shell_info_s close;
    close.shells_close_to_chunk = (int *)malloc(sizeof(int)*BasisSet.num_shells);
    close.aos_close_to_chunk = (int *)malloc(sizeof(int)*BasisSet.num_ao);
    close.close_shells_per_am = (int *)malloc(sizeof(int)*BasisSet.max_am);
    /*close.close_COCC = (double **)malloc(sizeof(double *)*BasisSet.num_ao);*/
    close.close_COCC = block_matrix(BasisSet.num_ao,MOInfo.ndocc);
    return close;
}

void free_close_shell_info(struct close_shell_info_s close_shell_info){
    
    free(close_shell_info.shells_close_to_chunk);
    free(close_shell_info.close_shells_per_am);
    free(close_shell_info.close_COCC);
    
}

void print_close_shell_info(struct close_shell_info_s close){
    
    int i;

    fprintf(outfile,"\nClose Shell Info Data Structure");
    
    fprintf(outfile,"\nNumber of close AO's = %d",close.num_close_aos);
    fprintf(outfile,"\nNumber of close Shells = %d",close.num_close_shells);
    for(i=0;i<BasisSet.num_shells;i++)
	fprintf(outfile,"\nshell_close_to_chunk[%d] = %d",i,close.shells_close_to_chunk[i]);
    for(i=0;i<BasisSet.num_ao;i++)
	fprintf(outfile,"\naos_close_to_chunk[%d] = %d",i,close.aos_close_to_chunk[i]);
    for(i=0;i<BasisSet.max_am;i++)
	fprintf(outfile,"\nclose_shells_per_am[%d] = %d",i,close.close_shells_per_am[i]);
    fprintf(outfile,"\nClose Occupied Eigenvector");
    print_mat(close.close_COCC,close.num_close_aos,MOInfo.ndocc,outfile);
    fprintf(outfile,"\n\n");
    fflush(outfile);
}
