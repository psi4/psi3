#define EXTERN
#include "includes.h"
#include "globals.h"
#include "prototypes.h"

void read_density()
{ 
  int i,j,k,l,dim_i,count,ibf;
  int *locs;
  PSI_FPTR ptr_start, junk;
  double **psq_so, *tmp_arr, **tmp_mat, **psq_ao;

  /* CDS 2/02 */
  int irrep, mo_offset, so_offset, i_ci, j_ci, max_opi, errcod;
  int fzc, populated_orbs, root;
  int *docc, *socc, *frozen_docc, *frozen_uocc, *reorder, **ras_opi; 
  double **scfvec, **opdm_blk, **onepdm;
  char *ref;

  Ptot = init_array(natri);

  rfile(opdm_file);

  /* Read in density in MO basis --- CDS 2/02 */
  if (!strcmp(opdm_basis,"MO")) {     

    errcod = ip_string("REFERENCE", &ref, 0);
    if (errcod == IPE_OK && !strcmp(ref, "UHF"))
      punt("Can't read MO OPDM for UHF yet");
          
    psq_so = block_matrix(nbfso,nbfso);

    max_opi = 0;
    for (irrep=0; irrep<nirreps; irrep++)
      if (sopi[irrep] > max_opi) max_opi = sopi[irrep];
    opdm_blk = block_matrix(max_opi, max_opi); 
    tmp_mat = block_matrix(max_opi, max_opi);

    docc = init_int_array(nirreps);
    socc = init_int_array(nirreps);
    frozen_docc = init_int_array(nirreps);
    frozen_uocc = init_int_array(nirreps);
    ras_opi = init_int_matrix(4,nirreps);
    reorder = init_int_array(nmo);

    fzc = 1;
    errcod = ip_boolean("FREEZE_CORE",&fzc,0);
 
    if (strcmp(wfn, "CI") == 0 || strcmp(wfn, "DETCI") == 0
      || strcmp(wfn, "DETCAS") == 0) {
      if (!ras_set(nirreps, nmo, fzc, orbspi, docc, socc,
		   frozen_docc, frozen_uocc, ras_opi, reorder, 1) )
        punt("Error in ras_set()");
    }
    else {
      errcod = ip_int_array("DOCC",docc,nirreps);
      if (errcod != IPE_OK) {
        free(docc);
	docc = chkpt_rd_clsdpi();
      }
      errcod = ip_int_array("SOCC",socc,nirreps);
      if (errcod != IPE_OK) {
        free(socc);
	socc = chkpt_rd_openpi();
      }

      reorder_qt(docc, socc, frozen_docc, frozen_uocc,
               reorder, orbspi, nirreps);
    }

    if (print_lvl >= PRINTOPDMLEVEL) {
      fprintf(outfile, "Reorder array:\n");
      for (i=0; i<nmo; i++) fprintf(outfile, "%d ", reorder[i]);
      fprintf(outfile, "\n");
    }

    populated_orbs = nmo;
    for (irrep=0; irrep<nirreps; irrep++)
      populated_orbs -= frozen_uocc[irrep];

    onepdm = block_matrix(populated_orbs, populated_orbs);
    root = 1;
    errcod = ip_data("ROOT","%d",&root,0);
    junk = 0;
    for (i=0; i<root; i++) {
      wreadw(opdm_file, (char *) onepdm[0], populated_orbs *
	     populated_orbs * sizeof(double), junk, &junk);
    }

    if (print_lvl > 2) {
      fprintf(outfile, "  Density matrix read in:\n");
      print_mat(onepdm,populated_orbs,populated_orbs,outfile);
      fprintf(outfile, "\n");
    }

    mo_offset = 0;
    so_offset = 0;
    for (irrep=0; irrep<nirreps; irrep++) {
      if (orbspi[irrep] == 0) continue;
      for (i=0; i<orbspi[irrep]-frozen_uocc[irrep]; i++) {
        for (j=0; j<orbspi[irrep]-frozen_uocc[irrep]; j++) {
          i_ci = reorder[i+mo_offset];
          j_ci = reorder[j+mo_offset];
          opdm_blk[i][j] = onepdm[i_ci][j_ci];
        }
      }

      if (print_lvl >= PRINTOPDMLEVEL) {
        fprintf(outfile, "Irrep %d (MO basis)\n", irrep);
        print_mat(opdm_blk,orbspi[irrep]-frozen_uocc[irrep], 
                  orbspi[irrep]-frozen_uocc[irrep],outfile);
      }

#if USE_LIBCHKPT
      scfvec = chkpt_rd_scf_irrep(irrep);
#else
      scfvec = file30_rd_blk_scf(irrep);
#endif
      mmult(opdm_blk,0,scfvec,1,tmp_mat,0,orbspi[irrep]-frozen_uocc[irrep],
            orbspi[irrep]-frozen_uocc[irrep],sopi[irrep],0);
      mmult(scfvec,0,tmp_mat,0,opdm_blk,0,sopi[irrep],
            orbspi[irrep]-frozen_uocc[irrep],sopi[irrep],0);
      if (print_lvl >= PRINTOPDMLEVEL) {
        fprintf(outfile,"Irrep %d (SO basis)\n", irrep);
        print_mat(opdm_blk,sopi[irrep], sopi[irrep],outfile);
      }
      for (i=0; i<sopi[irrep]; i++)
        for (j=0; j<sopi[irrep]; j++)
          psq_so[i+so_offset][j+so_offset] = opdm_blk[i][j];
      mo_offset += orbspi[irrep];
      so_offset += sopi[irrep];
    }

    if (print_lvl >= PRINTOPDMLEVEL) {
      fprintf(outfile,"  Total density matrix in SO basis :\n");
      print_mat(psq_so,nbfso,nbfso,outfile);
      fprintf(outfile,"\n");
    }
    
    free_block(onepdm);
    free_block(opdm_blk);
    free_block(tmp_mat);
    free_block(scfvec);
    free(docc); free(socc); free(frozen_docc); free(frozen_uocc);
    free(reorder);
    free_int_matrix(ras_opi, 4);

    tmp_mat = init_matrix(nbfso,nbfao);
    psq_ao = init_matrix(nbfao,nbfao);
    mmult(psq_so,0,usotao,0,tmp_mat,0,nbfso,nbfso,nbfao,0);
    mmult(usotao,1,tmp_mat,0,psq_ao,0,nbfao,nbfso,nbfao,0);
    free_matrix(tmp_mat,nbfso);
    free_block(psq_so);

    if (asymm_opdm && !strcmp(opdm_format,"SQUARE")) {
      for(i=0;i<nbfao;i++)
        for(j=0;j<=i;j++)
          Ptot[ioff[i]+j] = (psq_ao[i][j] + psq_ao[j][i])/2;
    }
    else
      sq_to_tri(psq_ao,Ptot,nbfao);

    free_matrix(psq_ao,nbfao);

  } /* end read MO opdm case */

  if (!strcmp(opdm_basis,"SO")) {	/* Read in density in SO basis */
    if (!strcmp(opdm_format,"SQUARE")) {
      psq_so = block_matrix(nbfso,nbfso);
      wreadw(opdm_file,(char *) psq_so[0],sizeof(double)*nbfso*nbfso,0,&junk);
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
      wreadw(opdm_file,(char *) psq_ao[0],sizeof(double)*nbfao*nbfao,0,&junk);
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

