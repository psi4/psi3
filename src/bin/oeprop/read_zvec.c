#define EXTERN
#include "includes.h"
#include "oeprop.h"
#include "oeprop.gbl"


void read_zvec()
{
  int i,j,m,n;
  PSI_FPTR ptr;
  double *tmp,sum;
  
  tmp = init_array(nstri);
  zvec = init_matrix(natri,natri);
  rfile(zvec_file);
  wreadw(zvec_file, (char *) tmp, sizeof(double)*nstri, 0, &ptr);
  if (delete_zvec)
    rclose(zvec_file,4);
  else
    rclose(zvec_file,3);
  if (print_lvl >= PRINTZVECTORLEVEL) {
    fprintf(outfile,"  Z-vector matrix in MO basis:\n");
    print_array(tmp,nbfso,outfile);
    fprintf(outfile,"\n\n");
  }

  for(m=0;m<nbfao;m++)
    for(n=0;n<nbfao;n++) {
      sum = 0;
      for(i=0;i<nbfso;i++)
        for(j=0;j<=i;j++)
          sum += scf_evec_ao[m][i]*scf_evec_ao[n][j]*tmp[ioff[i]+j];
      zvec[m][n] = sum;
    }
  
  free(tmp);

}
