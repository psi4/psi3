#define EXTERN
#define TOL 1E-14
/* #define DEBUG */

extern "C" {
   #include <stdio.h>
   #include <stdlib.h>
   /* may no longer need #include <libc.h> */
   #include <math.h>
   #include <libciomr.h>
   #include <qt.h>
   #include <file30.h>
   #include "structs.h"
   #include "globals.h"
}


#ifndef CIVECT_H
#include "civect.h"
#endif

void orbsfile_rd_blk(int targetfile, int root, int irrep, double **orbs_vector);
void orbsfile_wt_blk(int targetfile, int root, int irrep, double **orbs_vector);
void ave(int targetfile, double **onepdm);
void opdm_block(struct stringwr **alplist, struct stringwr **betlist,
		double **onepdm, double **CJ, double **CI, int Ja_list, 
		int Jb_list, int Jnas, int Jnbs, int Ia_list, int Ib_list, 
		int Inas, int Inbs);

void opdm(struct stringwr **alplist, struct stringwr **betlist, 
          int Inroots, int Iroot, int Inunits, int Ifirstunit, 
	  int Jnroots, int Jroot, int Jnunits, int Jfirstunit, 
	  int targetfile, int writeflag, int printflag)
{

  CIvect Ivec, Jvec;
  int i, j, k, l, maxrows, maxcols, bufsz, cino_index, roots;
  int max_orb_per_irrep;
  double **transp_tmp = NULL;
  double **transp_tmp2 = NULL;
  double *buffer1, *buffer2, **onepdm;
  int i_ci, j_ci, irrep, offset, orb_length=0, opdm_length=0;
  double *opdm_eigval, **opdm_eigvec, **opdm_blk, **scfvec, **scfvec30; 
  int Iblock, Iblock2, Ibuf, Iac, Ibc, Inas, Inbs, Iairr;
  int Jblock, Jblock2, Jbuf, Jac, Jbc, Jnas, Jnbs, Jairr;
  int do_Jblock, do_Jblock2;
  int populated_orbs;
  PSI_FPTR ljunk=0, onepdm_idx=0;
  double **tmp_mat, *tmp_mat2, **opdmso;

  Ivec.set(CIblks.vectlen, CIblks.num_blocks, Parameters.icore, Parameters.Ms0,
           CIblks.Ia_code, CIblks.Ib_code, CIblks.Ia_size, CIblks.Ib_size,
           CIblks.offset, CIblks.num_alp_codes, CIblks.num_bet_codes,
           CalcInfo.nirreps, AlphaG->subgr_per_irrep, Inroots, Inunits,
           Ifirstunit, CIblks.first_iablk, CIblks.last_iablk, CIblks.decode);

  Jvec.set(CIblks.vectlen, CIblks.num_blocks, Parameters.icore, Parameters.Ms0,
           CIblks.Ia_code, CIblks.Ib_code, CIblks.Ia_size, CIblks.Ib_size,
           CIblks.offset, CIblks.num_alp_codes, CIblks.num_bet_codes,
           CalcInfo.nirreps, AlphaG->subgr_per_irrep, Jnroots, Jnunits,
           Jfirstunit, CIblks.first_iablk, CIblks.last_iablk, CIblks.decode);

  if (Parameters.opdm_diag) rfile(Parameters.opdm_orbsfile);
  file30_init();

  populated_orbs = CalcInfo.num_ci_orbs + CalcInfo.num_fzc_orbs;
  for (irrep=0; irrep<CalcInfo.nirreps; irrep++) {
     opdm_length += (CalcInfo.orbs_per_irr[irrep] - CalcInfo.frozen_uocc[irrep])
                 * (CalcInfo.orbs_per_irr[irrep] - CalcInfo.frozen_uocc[irrep]);
     orb_length += (CalcInfo.orbs_per_irr[irrep]*CalcInfo.orbs_per_irr[irrep]);
     }
  max_orb_per_irrep = CalcInfo.max_orbs_per_irrep;
  opdm_eigvec = block_matrix(max_orb_per_irrep, max_orb_per_irrep);
  opdm_eigval = init_array(max_orb_per_irrep);
  opdm_blk = block_matrix(max_orb_per_irrep, max_orb_per_irrep);
  opdmso = block_matrix(CalcInfo.nbfso, CalcInfo.nbfso);
  tmp_mat = block_matrix(max_orb_per_irrep, max_orb_per_irrep);
  tmp_mat2 = init_array(2000);
  scfvec = block_matrix(max_orb_per_irrep, max_orb_per_irrep);
  Parameters.opdm_idxmat =
    init_int_matrix(Parameters.num_roots+2, CalcInfo.nirreps);
  Parameters.orbs_idxmat = 
    init_int_matrix(Parameters.num_roots+2, CalcInfo.nirreps);
  for (l=0; l<=(Parameters.num_roots+1); l++) {
     Parameters.opdm_idxmat[l][0] = l * opdm_length * sizeof(double); 
     Parameters.orbs_idxmat[l][0] = l * orb_length * sizeof(double);
     for (irrep=1; irrep<CalcInfo.nirreps; irrep++) {
        Parameters.orbs_idxmat[l][irrep] =
          Parameters.orbs_idxmat[l][irrep-1]+
          CalcInfo.orbs_per_irr[irrep-1]*CalcInfo.orbs_per_irr[irrep-1]
          *sizeof(double);  
        Parameters.opdm_idxmat[l][irrep] =
          Parameters.opdm_idxmat[l][irrep-1] +
          (CalcInfo.orbs_per_irr[irrep-1]-CalcInfo.frozen_uocc[irrep]) *
          (CalcInfo.orbs_per_irr[irrep-1]-CalcInfo.frozen_uocc[irrep]) *
          sizeof(double);
        } 
     } 

  buffer1 = Ivec.buf_malloc();
  buffer2 = Jvec.buf_malloc();
  Ivec.buf_lock(buffer1);
  Jvec.buf_lock(buffer2);
  onepdm = block_matrix(populated_orbs, populated_orbs); 

  if ((Ivec.icore==2 && Ivec.Ms0 && CalcInfo.ref_sym != 0) || 
      (Ivec.icore==0 && Ivec.Ms0)) {
    for (i=0, maxrows=0, maxcols=0; i<Ivec.num_blocks; i++) {
      if (Ivec.Ia_size[i] > maxrows) maxrows = Ivec.Ia_size[i];
      if (Ivec.Ib_size[i] > maxcols) maxcols = Ivec.Ib_size[i];
    }
    if (maxcols > maxrows) maxrows = maxcols;
    transp_tmp = (double **) malloc (maxrows * sizeof(double *));
    transp_tmp2 = (double **) malloc (maxrows * sizeof(double *));
    if (transp_tmp == NULL || transp_tmp2 == NULL) {
      printf("(opdm): Trouble with malloc'ing transp_tmp\n");
    }
    bufsz = Ivec.get_max_blk_size();
    transp_tmp[0] = init_array(bufsz);
    transp_tmp2[0] = init_array(bufsz);
    if (transp_tmp[0] == NULL || transp_tmp2[0] == NULL) {
      printf("(opdm): Trouble with malloc'ing transp_tmp[0]\n");
    }
  }
 
  for (k=0; k<Parameters.num_roots; k++) {
   
    if (k != 0) zero_mat(onepdm, populated_orbs, populated_orbs); 

    for (i=0; i<CalcInfo.num_fzc_orbs; i++)
      onepdm[i][i] = 2.0;

    if (Parameters.icore == 0) {
 
      for (Ibuf=0; Ibuf<Ivec.buf_per_vect; Ibuf++) {
        Ivec.read(Iroot, Ibuf);
        Iblock = Ivec.buf2blk[Ibuf];
        Iac = Ivec.Ia_code[Iblock];
        Ibc = Ivec.Ib_code[Iblock];
        Inas = Ivec.Ia_size[Iblock];
        Inbs = Ivec.Ib_size[Iblock];
       
        for (Jbuf=0; Jbuf<Jvec.buf_per_vect; Jbuf++) {
          do_Jblock=0; do_Jblock2=0;
          Jblock = Jvec.buf2blk[Jbuf];
          Jblock2 = -1;
          Jac = Jvec.Ia_code[Jblock];
          Jbc = Jvec.Ib_code[Jblock];
          if (Jvec.Ms0) Jblock2 = Jvec.decode[Jbc][Jac];
          Jnas = Jvec.Ia_size[Jblock];
          Jnbs = Jvec.Ib_size[Jblock];
          if (s1_contrib[Iblock][Jblock] || s2_contrib[Iblock][Jblock]) 
            do_Jblock = 1;
          if (Jvec.buf_offdiag[Jbuf] && (s1_contrib[Iblock][Jblock2] ||
                                         s2_contrib[Iblock][Jblock2]))
            do_Jblock2 = 1;
          if (!do_Jblock && !do_Jblock2) continue;
	 
          Jvec.read(Jroot, Jbuf);
	 
          if (do_Jblock) {
            opdm_block(alplist, betlist, onepdm, Jvec.blocks[Jblock], 
                       Ivec.blocks[Iblock], Jac, Jbc, Jnas,
                       Jnbs, Iac, Ibc, Inas, Inbs);
            }
	 
          if (do_Jblock2) {
            Jvec.transp_block(Jblock, transp_tmp);
            opdm_block(alplist, betlist, onepdm, transp_tmp,
                       Ivec.blocks[Iblock], Jbc, Jac, Jnbs,
                       Jnas, Iac, Ibc, Inas, Inbs);
          }
	 
       } /* end loop over Jbuf */
       
      if (Ivec.buf_offdiag[Ibuf]) { /* need to get contrib of transpose */
        Iblock2 = Ivec.decode[Ibc][Iac];
        Iac = Ivec.Ia_code[Iblock2];
        Ibc = Ivec.Ib_code[Iblock2];
        Inas = Ivec.Ia_size[Iblock2];
        Inbs = Ivec.Ib_size[Iblock2];
       
        Ivec.transp_block(Iblock, transp_tmp2);

        for (Jbuf=0; Jbuf<Jvec.buf_per_vect; Jbuf++) {
          do_Jblock=0; do_Jblock2=0;
          Jblock = Jvec.buf2blk[Jbuf];
          Jblock2 = -1;
          Jac = Jvec.Ia_code[Jblock];
          Jbc = Jvec.Ib_code[Jblock];
          if (Jvec.Ms0) Jblock2 = Jvec.decode[Jbc][Jac];
          Jnas = Jvec.Ia_size[Jblock];
          Jnbs = Jvec.Ib_size[Jblock];
          if (s1_contrib[Iblock2][Jblock] || s2_contrib[Iblock2][Jblock]) 
            do_Jblock = 1;
          if (Jvec.buf_offdiag[Jbuf] && (s1_contrib[Iblock2][Jblock2] ||
                                         s2_contrib[Iblock2][Jblock2]))
            do_Jblock2 = 1;
          if (!do_Jblock && !do_Jblock2) continue;
	   
          Jvec.read(Jroot, Jbuf);
	 
          if (do_Jblock) {
            opdm_block(alplist, betlist, onepdm, Jvec.blocks[Jblock], 
                       transp_tmp2, Jac, Jbc, Jnas,
                       Jnbs, Iac, Ibc, Inas, Inbs);
          }
	   
          if (do_Jblock2) {
            Jvec.transp_block(Jblock, transp_tmp);
            opdm_block(alplist, betlist, onepdm, transp_tmp,
                       transp_tmp2, Jbc, Jac, Jnbs,
                       Jnas, Iac, Ibc, Inas, Inbs);
          }
        } /* end loop over Jbuf */
      } /* end loop over Ibuf transpose */
    } /* end loop over Ibuf */
  } /* end icore==0 */

  else if (Parameters.icore==1) { /* whole vectors in-core */
    Ivec.read(Iroot, 0);
    Jvec.read(Jroot, 0);
    for (Iblock=0; Iblock<Ivec.num_blocks; Iblock++) {
      Iac = Ivec.Ia_code[Iblock];
      Ibc = Ivec.Ib_code[Iblock];
      Inas = Ivec.Ia_size[Iblock];
      Inbs = Ivec.Ib_size[Iblock];
      if (Inas==0 || Inbs==0) continue;
      for (Jblock=0; Jblock<Jvec.num_blocks; Jblock++) {
        Jac = Jvec.Ia_code[Jblock];
        Jbc = Jvec.Ib_code[Jblock];
        Jnas = Jvec.Ia_size[Jblock];
        Jnbs = Jvec.Ib_size[Jblock];
        if (s1_contrib[Iblock][Jblock] || s2_contrib[Iblock][Jblock])
          opdm_block(alplist, betlist, onepdm, Jvec.blocks[Jblock],
                     Ivec.blocks[Iblock], Jac, Jbc, Jnas,
                     Jnbs, Iac, Ibc, Inas, Inbs);
      }
    } /* end loop over Iblock */
  } /* end icore==1 */

  else if (Parameters.icore==2) { /* icore==2 */
    for (Ibuf=0; Ibuf<Ivec.buf_per_vect; Ibuf++) {
      Ivec.read(Iroot, Ibuf);
      Iairr = Ivec.buf2blk[Ibuf];

      for (Jbuf=0; Jbuf<Jvec.buf_per_vect; Jbuf++) {
        Jvec.read(Jroot, Jbuf);
        Jairr = Jvec.buf2blk[Jbuf];
	
      for (Iblock=Ivec.first_ablk[Iairr]; Iblock<=Ivec.last_ablk[Iairr];
           Iblock++) {
        Iac = Ivec.Ia_code[Iblock];
        Ibc = Ivec.Ib_code[Iblock];
        Inas = Ivec.Ia_size[Iblock];
        Inbs = Ivec.Ib_size[Iblock];
	   
        for (Jblock=Jvec.first_ablk[Jairr]; Jblock<=Jvec.last_ablk[Jairr];
             Jblock++) {
          Jac = Jvec.Ia_code[Jblock];
          Jbc = Jvec.Ib_code[Jblock];
          Jnas = Jvec.Ia_size[Jblock];
          Jnbs = Jvec.Ib_size[Jblock];
	   
          if (s1_contrib[Iblock][Jblock] || s2_contrib[Iblock][Jblock])
            opdm_block(alplist, betlist, onepdm, Jvec.blocks[Jblock],
                       Ivec.blocks[Iblock], Jac, Jbc, Jnas,
                       Jnbs, Iac, Ibc, Inas, Inbs);

          if (Jvec.buf_offdiag[Jbuf]) {
            Jblock2 = Jvec.decode[Jbc][Jac];
            if (s1_contrib[Iblock][Jblock2] ||
                s2_contrib[Iblock][Jblock2]) {
            Jvec.transp_block(Jblock, transp_tmp);
              opdm_block(alplist, betlist, onepdm, transp_tmp,
                Ivec.blocks[Iblock], Jbc, Jac,
                Jnbs, Jnas, Iac, Ibc, Inas, Inbs);
	    }
	  }

        } /* end loop over Jblock */

        if (Ivec.buf_offdiag[Ibuf]) {
          Iblock2 = Ivec.decode[Ibc][Iac];
          Ivec.transp_block(Iblock, transp_tmp2);
          Iac = Ivec.Ia_code[Iblock2];
          Ibc = Ivec.Ib_code[Iblock2];
          Inas = Ivec.Ia_size[Iblock2];
          Inbs = Ivec.Ib_size[Iblock2];
	   
          for (Jblock=Jvec.first_ablk[Jairr]; Jblock<=Jvec.last_ablk[Jairr];
            Jblock++) {
            Jac = Jvec.Ia_code[Jblock];
            Jbc = Jvec.Ib_code[Jblock];
            Jnas = Jvec.Ia_size[Jblock];
            Jnbs = Jvec.Ib_size[Jblock];
	   
            if (s1_contrib[Iblock2][Jblock] || s2_contrib[Iblock2][Jblock])
              opdm_block(alplist, betlist, onepdm, Jvec.blocks[Jblock],
                         transp_tmp2, Jac, Jbc, Jnas, Jnbs, Iac, Ibc, 
                         Inas, Inbs);

              if (Jvec.buf_offdiag[Jbuf]) {
                Jblock2 = Jvec.decode[Jbc][Jac];
                if (s1_contrib[Iblock][Jblock2] || 
                  s2_contrib[Iblock][Jblock2]) {
                  Jvec.transp_block(Jblock, transp_tmp);
                  opdm_block(alplist, betlist, onepdm, transp_tmp,
                    transp_tmp2, Jbc, Jac, Jnbs, Jnas, Iac, Ibc, Inas, Inbs);
                }
	      }

	    } /* end loop over Jblock */
          } /* end Ivec offdiag */

        } /* end loop over Iblock */
      } /* end loop over Jbuf */
    } /* end loop over Ibuf */
  } /* end icore==2 */

  else {
    printf("opdm: unrecognized core option!\n");
    return;
  }

  /* write and/or print the opdm */
  if (printflag) {
    fprintf(outfile, "\n\nOne-particle density matrix MO basis\n");
    print_mat(onepdm, populated_orbs, populated_orbs, outfile);
    fprintf(outfile, "\n");
  }


  if (writeflag) {
    rfile(targetfile);
    wwritw(targetfile, (char *) onepdm[0], populated_orbs *
      populated_orbs * sizeof(double), onepdm_idx, &onepdm_idx);
  }

   /* Convert the OPDMMO into pitzer ordering and backtransform to the SO */ 
   /* basis OPDMSO in pitzer ordering is written to the end of targetfile */
  /*
   offset = 0;
   fprintf(outfile,"OPDM in SO basis\n");
   for (irrep=0; irrep<CalcInfo.nirreps; irrep++) {
      if (CalcInfo.orbs_per_irr[irrep] == 0) continue;
      for (i=0; i<CalcInfo.orbs_per_irr[irrep]; i++) {
         for (j=0; j<CalcInfo.orbs_per_irr[irrep]; j++) {
            i_ci = CalcInfo.reorder[i+offset];
            j_ci = CalcInfo.reorder[j+offset];
            opdm_blk[i][j] = onepdm[i_ci][j_ci];
          }
        }
      scfvec30 = file30_rd_blk_scf(irrep);    
      mmult(opdm_blk,0,scfvec30,1,tmp_mat,0,CalcInfo.orbs_per_irr[irrep],
            CalcInfo.orbs_per_irr[irrep],CalcInfo.orbs_per_irr[irrep],0);
      mmult(scfvec30,0,tmp_mat,0,opdm_blk,0,CalcInfo.orbs_per_irr[irrep],
            CalcInfo.orbs_per_irr[irrep],CalcInfo.orbs_per_irr[irrep],0); 
      print_mat(opdm_blk,CalcInfo.orbs_per_irr[irrep],
                CalcInfo.orbs_per_irr[irrep],outfile);
      for (i=0; i<CalcInfo.orbs_per_irr[irrep]; i++) 
         for (j=0; j<CalcInfo.orbs_per_irr[irrep]; j++) 
            opdmso[i+offset][j+offset] = opdm_blk[i][j];
      offset += CalcInfo.orbs_per_irr[irrep];
     }

  if (writeflag) wwritw(targetfile, (char *) opdmso[0],sizeof(double)*
                         CalcInfo.nbfso*CalcInfo.nbfso, onepdm_idx, &onepdm_idx);
   */
  /*
  wreadw(targetfile, (char *) opdmso[0],sizeof(double)*
         CalcInfo.nbfso*CalcInfo.nbfso, 5000, &onepdm_idx);
   fprintf(outfile,"OPDM in SO basis read from file %d\n",targetfile);
     print_mat(opdmso,CalcInfo.nbfso,CalcInfo.nbfso,outfile);
  */

  if (writeflag) rclose(targetfile, 3);
  
  fflush(outfile);
  Iroot++; Jroot++;
  } /* end loop over num_roots k */  

  if (transp_tmp != NULL) free_block(transp_tmp);
  if (transp_tmp2 != NULL) free_block(transp_tmp2);
  Ivec.buf_unlock();
  Jvec.buf_unlock();
  free(buffer1);
  free(buffer2);


  /* Average the opdm's */
  if (Parameters.opdm_diag) rfile(targetfile);
  if (Parameters.opdm_ave) ave(Parameters.opdm_file, onepdm); 

  /* get CI Natural Orbitals */
  if (Parameters.opdm_diag) {

  /*   file30_init(); */

    /* reorder opdm from ci to pitzer and diagonalize each 
    ** symmetry blk
    */
    if (Parameters.opdm_ave) {
      cino_index = Parameters.num_roots;
      onepdm_idx = Parameters.num_roots*populated_orbs*
                   populated_orbs*sizeof(double);
      roots = Parameters.num_roots+1; 
    }
    else {
      cino_index = 0;
      onepdm_idx = 0;
      roots = Parameters.num_roots;
    }

    /* loop over roots or averaged opdm */
    for(k=cino_index; k<roots; k++) {
      if (Parameters.opdm_ave && Parameters.print_lvl > 1) {
        fprintf(outfile,"\n\n\t\t\tCI Natural Orbitals for the Averaged\n");
        fprintf(outfile,"\t\t\tOPDM of %d Roots in terms of Molecular"
                 " Orbitals\n\n",k); 
      }
      else if (Parameters.print_lvl > 1) {
        fprintf(outfile,
             "\n\t\t\tCI Natural Orbitals in terms of Molecular Orbitals\n\n");
        fprintf(outfile,"\t\t\t Root %d\n\n",k+1);
        fflush(outfile);
      }
      offset = 0;
      wreadw(targetfile, (char *) onepdm[0], populated_orbs * 
             populated_orbs * sizeof(double), onepdm_idx, &onepdm_idx); 
      for (irrep=0; irrep<CalcInfo.nirreps; irrep++) {
        if (CalcInfo.orbs_per_irr[irrep] == 0) continue; 
        for (i=0; i<CalcInfo.orbs_per_irr[irrep]-
                    CalcInfo.frozen_uocc[irrep]; i++) {
          for (j=0; j<CalcInfo.orbs_per_irr[irrep]-
                    CalcInfo.frozen_uocc[irrep]; j++) {
            i_ci = CalcInfo.reorder[i+offset];
            j_ci = CalcInfo.reorder[j+offset]; 
            opdm_blk[i][j] = onepdm[i_ci][j_ci];
          } 
        }
 
        if (k==0 || Parameters.opdm_ave) {
          /* Writting SCF vector to orbsfile for safe keeping */
          scfvec30 = file30_rd_blk_scf(irrep);
            #ifdef DEBUG
            fprintf(outfile,"Cvec for k==0, read in from file30 original\n");
            fprintf(outfile," %s Block \n", CalcInfo.labels[irrep]);
            print_mat(scfvec30, CalcInfo.orbs_per_irr[irrep],
                      CalcInfo.orbs_per_irr[irrep], outfile);
            #endif
          orbsfile_wt_blk(Parameters.opdm_orbsfile, Parameters.num_roots, 
                          irrep, scfvec30);
        }

        zero_mat(opdm_eigvec, max_orb_per_irrep, max_orb_per_irrep);

        if (CalcInfo.orbs_per_irr[irrep]-CalcInfo.frozen_uocc[irrep] > 0) {
          sq_rsp(CalcInfo.orbs_per_irr[irrep]-CalcInfo.frozen_uocc[irrep],
                 CalcInfo.orbs_per_irr[irrep]-CalcInfo.frozen_uocc[irrep],
                 opdm_blk, opdm_eigval, 1, opdm_eigvec, TOL); 
          }
        for (i=CalcInfo.orbs_per_irr[irrep]-CalcInfo.frozen_uocc[irrep]; 
             i<CalcInfo.orbs_per_irr[irrep]; i++) {
           opdm_eigvec[i][i] = 1.0;
           opdm_eigval[i] = 0.0;
           }
        eigsort(opdm_eigval, opdm_eigvec, -(CalcInfo.orbs_per_irr[irrep]));
        if (Parameters.print_lvl > 0) {
          if (irrep==0) {
            if (Parameters.opdm_ave) { 
              fprintf(outfile,
               "\n Averaged CI Natural Orbitals in terms of Molecular Orbitals\n\n");
              }
            else fprintf(outfile, 
                 "\n CI Natural Oribals in terms of Molecular Orbitals: Root %d\n\n",k+1);
          }
          fprintf(outfile,"\n %s Block \n", CalcInfo.labels[irrep]);        
          eivout(opdm_eigvec, opdm_eigval, CalcInfo.orbs_per_irr[irrep],
                 CalcInfo.orbs_per_irr[irrep], outfile);
        }

        if (Parameters.opdm_wrtnos) {
        /*
        */
          if (irrep==0) {
            if (Parameters.opdm_ave) {
              fprintf(outfile,"\n Writing CI Natural Orbitals from Averaged "
                "OPDM to orbsfile\n in terms of Symmetry Orbitals\n\n");
            }
            else {
              fprintf(outfile,"\n Writing CI Natural Orbitals for root %d"
                      " to orbsfile in terms of Symmetry Orbitals\n\n",k+1);
            }
          }
          orbsfile_rd_blk(Parameters.opdm_orbsfile, Parameters.num_roots, 
                          irrep, &scfvec[0]);  
          #ifdef DEBUG
          fprintf(outfile,"\nCvec read from orbsfile for MO to SO trans\n\n");
          fprintf(outfile," %s Block \n", CalcInfo.labels[irrep]);
          print_mat(scfvec, CalcInfo.orbs_per_irr[irrep],
                    CalcInfo.orbs_per_irr[irrep], outfile);
          fprintf(outfile,"\nOpdm_eigvec before MO to SO trans\n\n");
          fprintf(outfile," %s Block \n", CalcInfo.labels[irrep]);
          print_mat(opdm_eigvec, CalcInfo.orbs_per_irr[irrep],
                    CalcInfo.orbs_per_irr[irrep], outfile); 
          #endif
          mmult(scfvec, 0, opdm_eigvec, 0, opdm_blk, 0,
                CalcInfo.orbs_per_irr[irrep], CalcInfo.orbs_per_irr[irrep],
                CalcInfo.orbs_per_irr[irrep], 0); 
          if (Parameters.opdm_ave) { 
            orbsfile_wt_blk(Parameters.opdm_orbsfile, Parameters.num_roots+1, 
                            irrep, opdm_blk); 
          }
          else orbsfile_wt_blk(Parameters.opdm_orbsfile, k, irrep, opdm_blk); 
          #ifdef DEBUG
          fprintf(outfile,"\nOpdm_blk right after writting to orbsfile\n\n");
          fprintf(outfile," %s Block \n", CalcInfo.labels[irrep]);
          print_mat(opdm_blk, CalcInfo.orbs_per_irr[irrep],
                    CalcInfo.orbs_per_irr[irrep], outfile); 
          #endif
        }
        offset += CalcInfo.orbs_per_irr[irrep];
      } /* end loop over irreps */
    } /* end loop over roots */

    free(scfvec30);
    free_block(onepdm);

    /* now get appropriate orbitals from orbsfile and write to file30 */

    if (Parameters.opdm_wrtnos) {
      if (Parameters.opdm_ave) {
        for (irrep=0; irrep<CalcInfo.nirreps; irrep++) {
          if (irrep==0) {
            fprintf(outfile,"\n Writing CI Natural Orbitals from Averaged "
                "OPDM to orbsfile\n in terms of Symmetry Orbitals\n\n");
            }
          fprintf(outfile," %s Block \n", CalcInfo.labels[irrep]);
          orbsfile_rd_blk(Parameters.opdm_orbsfile, Parameters.num_roots+1, 
                          irrep, opdm_blk); 
          print_mat(opdm_blk, CalcInfo.orbs_per_irr[irrep],
                    CalcInfo.orbs_per_irr[irrep], outfile);
          file30_wt_blk_scf(opdm_blk, irrep);
          fprintf(outfile, "\n Warning: Natural Orbitals for the Averaged ");
          fprintf(outfile, "OPDM Have Been Written to file30!\n\n"); 
        }
      }
      else {
        for(k=0; k<Parameters.num_roots; k++) {
          for(irrep=0; irrep<CalcInfo.nirreps; irrep++) {
            if (irrep==0) {
              fprintf(outfile,"\n Writing CI Natural Orbitals for root %d"
                      " to orbsfile in terms of Symmetry Orbitals\n\n",k+1);
              }
            fprintf(outfile,"\n %s Block \n", CalcInfo.labels[irrep]);
            orbsfile_rd_blk(Parameters.opdm_orbsfile, k, irrep, opdm_blk); 
            print_mat(opdm_blk, CalcInfo.orbs_per_irr[irrep],
                      CalcInfo.orbs_per_irr[irrep], outfile);
            if (k==Parameters.opdm_orbs_root) { 
              file30_wt_blk_scf(opdm_blk, irrep);
              fprintf(outfile,"\n Warning: Natural Orbitals Have Been "
                    "Written to file30!\n\n"); 
            }
          } /* end irrep */
        } /* end num_roots */
      } /* end else */
    }
  } 
  /* CINOS completed */
 
  file30_close();

  fflush(outfile);
  if (Parameters.opdm_diag) {
    rclose(targetfile, 3);
    rclose(Parameters.opdm_orbsfile, 3);
    free_block(scfvec);
    free_block(opdm_blk);
    free_block(opdm_eigvec);
    free(opdm_eigval); 
  }

}

void opdm_block(struct stringwr **alplist, struct stringwr **betlist,
		double **onepdm, double **CJ, double **CI, int Ja_list, 
		int Jb_list, int Jnas, int Jnbs, int Ia_list, int Ib_list, 
		int Inas, int Inbs)
{
  int Ia_idx, Ib_idx, Ja_idx, Jb_idx, Ja_ex, Jb_ex, Jbcnt, Jacnt; 
  struct stringwr *Jb, *Ja;
  signed char *Jbsgn, *Jasgn;
  unsigned int *Jbridx, *Jaridx;
  double C1, C2, Ib_sgn, Ia_sgn;
  int i, j, oij, nfzc, *Jboij, *Jaoij;
 
  nfzc = CalcInfo.num_fzc_orbs;

  /* loop over Ia in Ia_list */
  if (Ia_list == Ja_list) {
    for (Ia_idx=0; Ia_idx<Inas; Ia_idx++) {
      for (Jb=betlist[Jb_list], Jb_idx=0; Jb_idx<Jnbs; Jb_idx++, Jb++) {
	C1 = CJ[Ia_idx][Jb_idx];

	/* loop over excitations E^b_{ij} from |B(J_b)> */
	Jbcnt = Jb->cnt[Ib_list];
	Jbridx = Jb->ridx[Ib_list];
	Jbsgn = Jb->sgn[Ib_list];
	Jboij = Jb->oij[Ib_list];
	for (Jb_ex=0; Jb_ex < Jbcnt; Jb_ex++) {
	  oij = *Jboij++;
	  Ib_idx = *Jbridx++;
	  Ib_sgn = (double) *Jbsgn++;
	  C2 = CI[Ia_idx][Ib_idx];
          i = oij/CalcInfo.num_ci_orbs + nfzc;
          j = oij%CalcInfo.num_ci_orbs + nfzc;
	  onepdm[i][j] += C1 * C2 * Ib_sgn;
	}
      }
    }
  }

  /* loop over Ib in Ib_list */
  if (Ib_list == Jb_list) {
    for (Ib_idx=0; Ib_idx<Inbs; Ib_idx++) {
      for (Ja=alplist[Ja_list], Ja_idx=0; Ja_idx<Jnas; Ja_idx++, Ja++) {
	C1 = CJ[Ja_idx][Ib_idx];
	
	/* loop over excitations */
	Jacnt = Ja->cnt[Ia_list];
	Jaridx = Ja->ridx[Ia_list];
	Jasgn = Ja->sgn[Ia_list];
	Jaoij = Ja->oij[Ia_list];
	for (Ja_ex=0; Ja_ex < Jacnt; Ja_ex++) {
	  oij = *Jaoij++;
	  Ia_idx = *Jaridx++;
	  Ia_sgn = (double) *Jasgn++;
	  C2 = CI[Ia_idx][Ib_idx];
          i = oij/CalcInfo.num_ci_orbs + nfzc;
          j = oij%CalcInfo.num_ci_orbs + nfzc;
	  onepdm[i][j] += C1 * C2 * Ia_sgn;
	}
      }
    }
  }
}



/*
** orbsfile_wt_blk()
**
** Parameters:
**    root        = number of the current CI root 
**    irrep       = number of the current irrep
**    orbs_vector = orbs_vector containing the block of orbitals
**                  for the current root and irrep number 
*/
void orbsfile_wt_blk(int targetfile, int root, int irrep, double **orbs_vector)
{
  int i, j, count;
  double *tmp_vector;
  PSI_FPTR jnk=0;

  if (CalcInfo.orbs_per_irr[irrep])
    {
    tmp_vector = init_array(CalcInfo.orbs_per_irr[irrep]*
                            CalcInfo.orbs_per_irr[irrep]);
    count = 0;
    for(i=0; i<CalcInfo.orbs_per_irr[irrep]; i++)
       for(j=0; j<CalcInfo.orbs_per_irr[irrep]; j++, count++) {
          tmp_vector[count] = orbs_vector[j][i];
         }
    wwritw(targetfile, (char *) tmp_vector, sizeof(double)*
           CalcInfo.orbs_per_irr[irrep]*CalcInfo.orbs_per_irr[irrep],
           (PSI_FPTR) Parameters.orbs_idxmat[root][irrep], &jnk);
    free(tmp_vector);
    }
}



/*
** orbsfile_rd_blk()
**
** Parameters:
**    root   = number of the current CI root
**    irrep  = number of the current irrep
*/
void orbsfile_rd_blk(int targetfile, int root, int irrep, double **orbs_vector)
{
  int i, j, count;
  double *tmp_vector;
  PSI_FPTR jnk=0;

  if(CalcInfo.orbs_per_irr[irrep])
    {
    tmp_vector = init_array(CalcInfo.orbs_per_irr[irrep]*
                            CalcInfo.orbs_per_irr[irrep]);

    wreadw(targetfile, (char *) tmp_vector, sizeof(double)*
           CalcInfo.orbs_per_irr[irrep]*CalcInfo.orbs_per_irr[irrep],
           (PSI_FPTR) Parameters.orbs_idxmat[root][irrep], &jnk);

    count = 0;
    for(i=0; i<CalcInfo.orbs_per_irr[irrep]; i++)
      for(j=0; j<CalcInfo.orbs_per_irr[irrep]; j++, count++) {
         orbs_vector[j][i] = tmp_vector[count];
      } 
    free(tmp_vector);
    }
}


/*
** ave()
** 
** Parameters:
**  targetfile = file number to obtain matrices from 
**  opdm_blk = use old opdm_blk matrix to store tmp matrices
**  opdm_eigvec = use old opdm_eigvec matrix to store tmp matrices 
**  max_orb_per_irrep = maximum orbital for all irreps
**
*/
void ave(int targetfile, double **tmp_mat1)
{
  int root, i, j, populated_orbs;
  double **tmp_mat2;
  PSI_FPTR jnk;

  populated_orbs = CalcInfo.nbfso-CalcInfo.num_fzv_orbs;
  zero_mat(tmp_mat1, populated_orbs, populated_orbs);
  tmp_mat2 = block_matrix(populated_orbs, populated_orbs);  

  for(root=0; root<Parameters.num_roots; root++) {

    if (root==0) {
      wreadw(targetfile,(char *) tmp_mat1[0], 
             populated_orbs * populated_orbs * sizeof(double), 
             0, &jnk);
      if (Parameters.opdm_print) {
        fprintf(outfile,"\n\n\t\tOPDM for Root 1");
        print_mat(tmp_mat1, populated_orbs, populated_orbs, outfile);
      }
    }

    else {
      wreadw(targetfile, (char *) tmp_mat2[0], 
             populated_orbs*populated_orbs*sizeof(double),
             root*populated_orbs*populated_orbs*sizeof(double), &jnk);
      for(i=0; i<populated_orbs; i++)
         for(j=0; j<populated_orbs; j++) 
            tmp_mat1[i][j] += tmp_mat2[i][j];    
      if (Parameters.opdm_print) {
        fprintf(outfile,"\n\n\t\tOPDM for Root %d",root+1);
        print_mat(tmp_mat2, populated_orbs, populated_orbs, outfile);
      }
    }

    zero_mat(tmp_mat2, populated_orbs, populated_orbs);
  }

  free(tmp_mat2);
  for (i=0; i<populated_orbs; i++)
    for (j=0; j<populated_orbs; j++) 
       tmp_mat1[i][j] *= (1.0/((double)Parameters.num_roots));
       wwritw(targetfile, (char *) tmp_mat1[0], 
        populated_orbs*populated_orbs*sizeof(double), 
         (Parameters.num_roots)*populated_orbs*populated_orbs*sizeof(double),&jnk);
  if (Parameters.print_lvl > 0 || Parameters.opdm_print) {
    fprintf(outfile,
      "\n\t\t\t Averaged OPDM's for %d Roots written to opdm_file \n\n",
      Parameters.num_roots);
  }
  if (Parameters.opdm_print) {
    print_mat(tmp_mat1, populated_orbs, populated_orbs, outfile);
  }

}


