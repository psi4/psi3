#define EXTERN
#include "includes.h"
#include "oeprop.gbl"
#include "oeprop.h"

void read_density()
{ 
  int i,j,k,l,dim_i,count,ibf;
  int *locs;
  PSI_FPTR ptr_start, junk;
  double **psq_so, *tmp_arr, **tmp_mat, **psq_ao;

  Ptot = init_array(natri);

  rfile(opdm_file);

  if (!strcmp(opdm_basis,"SO")) {	/* Read in density in SO basis */
    if (!strcmp(opdm_format,"SQUARE")) {
      psq_so = block_matrix(nbfso,nbfso);
      wreadw(opdm_file,(char *) psq_so[0], sizeof(double)*nbfso*nbfso, 0, &junk);
    }
    else {
      tmp_arr = init_array(nstri);
      psq_so = init_matrix(nbfso,nbfso);
      wreadw(opdm_file,(char *) tmp_arr, sizeof(double)*nstri, 0, &junk);
      tri_to_sq(tmp_arr,psq_so,nbfso);
      free(tmp_arr);
    }
          
          
    tmp_mat = init_matrix(nbfso,nbfao);
    psq_ao = init_matrix(nbfao,nbfao);
    mmult(psq_so,0,usotao,0,tmp_mat,0,nbfso,nbfso,nbfao,0);
    mmult(usotao,1,tmp_mat,0,psq_ao,0,nbfao,nbfso,nbfao,0);
    free_matrix(tmp_mat,nbfso);

    if (asymm_opdm && !strcmp(opdm_format,"SQUARE")) {
      for(i=0;i<nbfao;i++)
        for(j=0;j<=i;j++)
          Ptot[ioff[i]+j] = (psq_ao[i][j] + psq_ao[j][i])/2;
    }
    else
      sq_to_tri(psq_ao,Ptot,nbfao);

    free_matrix(psq_ao,nbfao);
    if (print_lvl >= PRINTOPDMLEVEL) {
      fprintf(outfile,"  Total density matrix in SO basis :\n");
      print_mat(psq_so,nbfso,nbfso,outfile);
      fprintf(outfile,"\n");
    }
    free_block(psq_so);
  }
  else if (!strcmp(opdm_basis,"AO")) {	/* Read in density in AO basis */
    if (!strcmp(opdm_format,"SQUARE")) {
      psq_ao = block_matrix(nbfao,nbfao);
      wreadw(opdm_file,(char *) psq_ao[0], sizeof(double)*nbfao*nbfao, 0, &junk);
      if (asymm_opdm) {
        for(i=0;i<nbfao;i++)
          for(j=0;j<=i;j++)
            Ptot[ioff[i]+j] = (psq_ao[i][j] + psq_ao[j][i])/2;
      }
      else
        sq_to_tri(psq_ao,Ptot,nbfao);
      free_matrix(psq_ao,nbfao);
    }
    else
      wreadw(opdm_file,(char *) Ptot, sizeof(double)*natri, 0, &junk);
  }

  rclose(opdm_file,3);
  
  if (print_lvl >= PRINTOPDMLEVEL) {
    fprintf(outfile,"  Total density matrix in AO basis :\n");
    print_array(Ptot,nbfao,outfile);
    fprintf(outfile,"\n\n");
  }

  return;
}

