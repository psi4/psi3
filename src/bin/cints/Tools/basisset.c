#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<libciomr.h>
#include<file30.h>

#include<libint.h>
#include"defines.h"
#define EXTERN
#include"global.h"
#include"shell_pairs.h"
#include"small_fns.h"

/*--------------------------------
  Explicit functions declarations
 --------------------------------*/
static void get_primitives(void);
static void get_shell_info(void);

void init_basisset()
{
  BasisSet.num_shells = file30_rd_nshell();
  BasisSet.num_prims = file30_rd_nprim();
  BasisSet.num_ao = file30_rd_nao();
  BasisSet.puream = (Symmetry.num_so != BasisSet.num_ao) ? 1 : 0;  /* need to transform to pure. ang. mom. basis? */
/* BasisSet.cgtos = */ get_primitives();
/* BasisSet.shells = */ get_shell_info();
/* BasisSet.shell_pairs = */ init_shell_pairs();

  /*-----------------------------------------------
    Namespaces are not well defined in CINTS,
    because of this here're some overlapping inits
   -----------------------------------------------*/
  if (BasisSet.puream && UserOptions.symm_ints)
    Symmetry.usotao = file30_rd_usotbf();
  else
    Symmetry.usotao = file30_rd_usotao_new();
  if (Symmetry.nirreps > 1 && UserOptions.symm_ints)
/* Symmetry.us_pairs = */ init_unique_shell_pairs();

  return;
}


void cleanup_basisset()
{
  int i;
  
  dealloc_pairs();
  for(i=0;i<BasisSet.num_shells;i++)
    free(BasisSet.shells[i].trans_vec);
  free(BasisSet.shells);
  free(BasisSet.cgtos);

  return;
}

void get_shell_info()
{
   int i, j, l, g, count, stab_index;
   int *shell_center;		/* atomic center of each shell */
   int *shell_type;		/* angular mom. of shell */
   int *shell_num_prims;	/* number of primitives per shell */
   int *prim_pointers;		/* first primitive in shell */
   int *shell_fbf;              /* first basisfn in shell */
   int *shell_fao;              /* first AO in shell */
   int **shell_trans_table;     /* shell transformation table */

   /*--- retrieve location of shells (which atom it's centered on) ---*/
   shell_center = file30_rd_snuc();

   /*--- retrieve angular momentum of each shell (1=s, 2=p, 3=d, etc  ) ---*/
   shell_type = file30_rd_stype();

   /*--- retrieve number of primitives per shell ---*/
   shell_num_prims = file30_rd_snumg();

   /*--- retrieve pointer to first primitive in shell ---*/
   prim_pointers = file30_rd_sprim();

   /*--- retrieve pointer to first basisfn in shell ---*/
   shell_fbf = file30_rd_sloc_new();

   /*--- retrieve pointer to first AO in shell ---*/
   shell_fao = file30_rd_sloc();
   
   /*--- retrieve shell tranformation table ---*/
   shell_trans_table = file30_rd_shell_transm();

   
   BasisSet.shells = (struct shell_def *) malloc(sizeof(struct shell_def)*
						 BasisSet.num_shells);
   BasisSet.max_am = 0;
   BasisSet.max_num_prims = 0;
   for (i=0; i<BasisSet.num_shells; i++){
      BasisSet.shells[i].center = shell_center[i];
      if ( (l = shell_type[i]) <= MAX_AM)
	BasisSet.shells[i].am = shell_type[i];
      else
	punt("Angular momentum limit exceeded");
      if (l > BasisSet.max_am)
        BasisSet.max_am = l;
      BasisSet.shells[i].n_prims = shell_num_prims[i];
      if (shell_num_prims[i] > BasisSet.max_num_prims)
        BasisSet.max_num_prims = shell_num_prims[i];
      BasisSet.shells[i].fprim = prim_pointers[i];
      BasisSet.shells[i].trans_vec = shell_trans_table[i];
      BasisSet.shells[i].fbf = shell_fbf[i];
      BasisSet.shells[i].fao = shell_fao[i];
      /*--- compute index of the stabilizer for the shell ---*/
      count = 1;
      for(g=1;g<Symmetry.nirreps;g++)
	if (i == BasisSet.shells[i].trans_vec[g]-1)
	  count++;
      stab_index = Symmetry.nirreps/count;
      if (Symmetry.max_stab_index < stab_index)
	Symmetry.max_stab_index = stab_index;
   }

   free(shell_center);
   free(shell_type);
   free(shell_num_prims);
   free(prim_pointers);
   free(shell_fbf);
   free(shell_fao);
   
   return;
}


void get_primitives(void)
{
   int i, j;
   double *exponents;      /* primitive gaussian exponents */
   double **ccoeffs;       /* primitive gaussian cont. coeff. for each ang. momentum*/

   /*--- read in exponents of primitive gaussians ---*/
   exponents = file30_rd_exps();

   /*--- read in coefficients of primitive gaussians ---*/
   ccoeffs = file30_rd_contr_full();

   /*--- allocate prims structure ---*/
   BasisSet.cgtos = (struct gaussian_function *)malloc(sizeof(struct gaussian_function)*BasisSet.num_prims);

   /*--- fill prims structure ---*/
   for (i=0; i<BasisSet.num_prims; i++){
     BasisSet.cgtos[i].exp = exponents[i];
     for(j=0;j<MAX_AM;j++) 
       BasisSet.cgtos[i].ccoeff[j] = ccoeffs[i][j];
   }

   free(exponents);
   free_matrix(ccoeffs,BasisSet.num_prims);

   return;
}


