#include <math.h>
#include <stdio.h>
#include <ip_libv1.h>
#include <libciomr.h>
#include <libint.h>
#include "build_libderiv.h"

FILE *infile, *outfile, *d1hrr_header,
     *deriv_header, *libderiv_header, *init_code;
int libderiv_stack_size[MAX_AM/2+1];
LibderivParams_t Params;

void punt();
extern void emit_deriv1_managers();
extern void emit_deriv2_managers();
extern void emit_deriv12_managers();
extern int emit_d1hrr_build();
extern int emit_d1hrr_build_macro();
extern int emit_deriv_build();
extern int emit_deriv_build_macro();

int main()
{
  int i,j,k,l,f;
  int j_min, j_max, k_min, k_max, l_min, l_max;
  int errcod;
  int new_am;
  int class_size;
  int num_subfunctions;
  int max_class_size = 785;
  const int io[] = {0,1,3,6,10,15,21,28,36,45,55,66,78,91,105,120,136,153};
  const char am_letter[] = "0pdfghiklmnoqrtuvwxyz";


  /*-------------------------------
    Initialize files and libraries
   -------------------------------*/
  infile = fopen("./input.dat", "r");
  outfile = fopen("./output.dat", "w");
  d1hrr_header = fopen("./d1hrr_header.h","w");
  deriv_header = fopen("./deriv_header.h","w");
  libderiv_header = fopen("./libderiv.h","w");
  init_code = fopen("./init_libderiv.cc","w");

  ip_set_uppercase(1);
  ip_initialize(infile,outfile);
  ip_cwk_add(":DEFAULT");
  ip_cwk_add(":LIBDERIV");
  
  /*---------------------------------------------
    Getting the new_am and deriv_lvl from user
    and making sure it is consistent with libint.h
   ---------------------------------------------*/
  errcod = ip_data("NEW_AM","%d",&new_am,0);
  if (errcod != IPE_OK)
    new_am = DEFAULT_NEW_AM;
  if (new_am <= 0)
    punt("  NEW_AM must be positive.");
  if (new_am > (LIBINT_MAX_AM - 1)*2-DERIV_LVL)
    punt("  Maximum NEW_AM is greater than the installed libint.a allows.\n  Recompile libint.a with greater NEW_AM.");

  /*-------------
    Init globals
   -------------*/
  for(l=0;l<=new_am/2;l++)
    libderiv_stack_size[l] = 0;
  Params.new_am = new_am;
  Params.old_am = 0;
  Params.opt_am = LIBINT_OPT_AM;
  Params.max_am_to_inline_vrr_worker = -1;
  Params.max_am_manager_to_inline_vrr_worker = -1;
  Params.max_am_to_inline_deriv_worker = -1;
  Params.max_am_manager_to_inline_deriv_worker = -1;
  Params.max_am_to_inline_hrr_worker = -1;
  Params.max_am_manager_to_inline_hrr_worker = -1;
  Params.max_am_to_inline_d1hrr_worker = -1;
  Params.max_am_manager_to_inline_d1hrr_worker = -1;
  Params.max_am_to_inline_vrr_manager = -1;

  /* Setting up init_libderiv.c, header.h */
  fprintf(init_code,"#include <stdlib.h>\n");
  fprintf(init_code,"#include <strings.h>\n");
  fprintf(init_code,"#include <libint.h>\n");
  fprintf(init_code,"#include \"libderiv.h\"\n");
  fprintf(init_code,"#include \"d1hrr_header.h\"\n\n");
  fprintf(init_code,"extern \"C\" {\n");
  fprintf(init_code,"/* This function initializes a matrix of pointers to routines */\n");
  fprintf(init_code,"/* for computing first derivatives of ERI classes up to (%cs|%cs) */\n",am_letter[new_am],am_letter[new_am]);
  fprintf(init_code,"void (*build_deriv1_eri[%d][%d][%d][%d])(Libderiv_t *, int);\n\n",new_am/2+1,new_am/2+1,new_am/2+1,new_am/2+1);
  fprintf(init_code,"/* This function initializes a matrix of pointers to routines */\n");
  fprintf(init_code,"/* for computing second derivatives of ERI classes up to (%cs|%cs) */\n",am_letter[new_am],am_letter[new_am]);
  fprintf(init_code,"void (*build_deriv2_eri[%d][%d][%d][%d])(Libderiv_t *, int);\n\n",new_am/2+1,new_am/2+1,new_am/2+1,new_am/2+1);
  fprintf(init_code,"/* This function initializes a matrix of pointers to routines */\n");
  fprintf(init_code,"/* for computing first and second derivatives of ERI classes up to (%cs|%cs) */\n",am_letter[new_am],am_letter[new_am]);
  fprintf(init_code,"void (*build_deriv12_eri[%d][%d][%d][%d])(Libderiv_t *, int);\n\n",new_am/2+1,new_am/2+1,new_am/2+1,new_am/2+1);
  fprintf(init_code,"int libderiv_stack_size[%d];\n",new_am/2+1);
  fprintf(init_code,"void init_libderiv_base()\n{\n");

  emit_deriv1_managers();
  emit_deriv12_managers();
#if EMIT_DERIV2_MANAGERS
  emit_deriv2_managers();
#endif
  emit_d1hrr_build();
  emit_d1hrr_build_macro();
  emit_deriv_build();
  emit_deriv_build_macro();

  /* put computed stack sizes for each angular momentum level into init_libderiv_base() */
  for(l=0;l<=new_am/2;l++)
    fprintf(init_code,"\n  libderiv_stack_size[%d] = %d;",l,libderiv_stack_size[l]);
  
  fprintf(init_code,"\n}\n\n");
  fprintf(init_code,"/* These functions initialize library objects */\n");
  fprintf(init_code,"/* Library objects operate independently of each other */\n");
  fprintf(init_code,"int init_libderiv(Libderiv_t *libderiv, int max_am, int max_num_prim_quartets, int max_cart_class_size)\n{\n");
  fprintf(init_code,"  int memory = 0;\n\n");
  fprintf(init_code,"  if (max_am >= LIBDERIV_MAX_AM) return -1;\n");
  fprintf(init_code,"  libderiv->int_stack = (double *) malloc(libderiv_stack_size[max_am]*sizeof(double));\n");
  fprintf(init_code,"  memory += libderiv_stack_size[max_am];\n");
  fprintf(init_code,"  libderiv->zero_stack = (double *) malloc(max_cart_class_size*sizeof(double));\n");
  fprintf(init_code,"  bzero((char *)libderiv->zero_stack,max_cart_class_size*sizeof(double));\n");
  fprintf(init_code,"  memory += max_cart_class_size;\n");
  fprintf(init_code,"  libderiv->PrimQuartet = (prim_data *) malloc(max_num_prim_quartets*sizeof(prim_data));\n");
  fprintf(init_code,"  memory += max_num_prim_quartets*sizeof(prim_data)/sizeof(double);\n");
  fprintf(init_code,"  return memory;\n}\n\n");
  fprintf(init_code,"void free_libderiv(Libderiv_t *libderiv)\n{\n");
  fprintf(init_code,"  if (libderiv->int_stack != NULL) {\n");
  fprintf(init_code,"    free(libderiv->int_stack);\n");
  fprintf(init_code,"    libderiv->int_stack = NULL;\n");
  fprintf(init_code,"  }\n");
  fprintf(init_code,"  if (libderiv->zero_stack != NULL) {\n");
  fprintf(init_code,"    free(libderiv->zero_stack);\n");
  fprintf(init_code,"    libderiv->zero_stack = NULL;\n");
  fprintf(init_code,"  }\n");
  fprintf(init_code,"  if (libderiv->PrimQuartet != NULL) {\n");
  fprintf(init_code,"    free(libderiv->PrimQuartet);\n");
  fprintf(init_code,"    libderiv->PrimQuartet = NULL;\n");
  fprintf(init_code,"  }\n\n");
  fprintf(init_code,"  return;\n}\n\n");
  fprintf(init_code,"int libderiv_storage_required(int max_am, int max_num_prim_quartets, int max_cart_class_size)\n{\n");
  fprintf(init_code,"  int memory = 0;\n\n");
  fprintf(init_code,"  if (max_am >= LIBDERIV_MAX_AM) return -1;\n");
  fprintf(init_code,"  memory += libderiv_stack_size[max_am];\n");
  fprintf(init_code,"  memory += max_cart_class_size;\n");
  fprintf(init_code,"  memory += max_num_prim_quartets*sizeof(prim_data)/sizeof(double);\n");
  fprintf(init_code,"  return memory;\n}\n");
  fprintf(init_code,"}\n"); /* end of extern "C" */
  fclose(init_code);
  fclose(d1hrr_header);
  fclose(deriv_header);
  
    /* Setting up libderiv.h */
  fprintf(libderiv_header,"#ifndef _psi3_libderiv_h\n");
  fprintf(libderiv_header,"#define _psi3_libderiv_h\n\n");
  fprintf(libderiv_header,"#include <libint.h>\n\n");
  fprintf(libderiv_header,"/* Maximum angular momentum of functions in a basis set plus 1 */\n");
  fprintf(libderiv_header,"#define LIBDERIV_MAX_AM %d\n",1+new_am/2);
  fprintf(libderiv_header,"#ifdef DERIV_LVL\n");
  fprintf(libderiv_header," #undef DERIV_LVL\n");
  fprintf(libderiv_header,"#endif\n");
  fprintf(libderiv_header,"#define DERIV_LVL %d\n\n",DERIV_LVL);
  fprintf(libderiv_header,"typedef struct {\n");
  fprintf(libderiv_header,"  double *int_stack;\n"); 
  fprintf(libderiv_header,"  prim_data *PrimQuartet;\n");
  fprintf(libderiv_header,"  double *zero_stack;\n");
  fprintf(libderiv_header,"  double *ABCD[12+144];\n");
  fprintf(libderiv_header,"  double AB[3];\n");
  fprintf(libderiv_header,"  double CD[3];\n");
  fprintf(libderiv_header,"  double *deriv_classes[%d][%d][%d];\n",1+new_am,1+new_am,12);
  fprintf(libderiv_header,"  double *deriv2_classes[%d][%d][%d];\n",1+new_am,1+new_am,144);
  fprintf(libderiv_header,"  double *dvrr_classes[%d][%d];\n",1+new_am,1+new_am);
  fprintf(libderiv_header,"  double *dvrr_stack;\n");
  fprintf(libderiv_header,"  } Libderiv_t;\n\n");
  fprintf(libderiv_header,"#ifdef __cplusplus\n");
  fprintf(libderiv_header,"extern \"C\" {\n");
  fprintf(libderiv_header,"#endif\n");
  fprintf(libderiv_header,"extern void (*build_deriv1_eri[%d][%d][%d][%d])(Libderiv_t *, int);\n",
	  new_am/2+1,new_am/2+1,new_am/2+1,new_am/2+1);
  fprintf(libderiv_header,"extern void (*build_deriv2_eri[%d][%d][%d][%d])(Libderiv_t *, int);\n",
	  new_am/2+1,new_am/2+1,new_am/2+1,new_am/2+1);
  fprintf(libderiv_header,"extern void (*build_deriv12_eri[%d][%d][%d][%d])(Libderiv_t *, int);\n",
	  new_am/2+1,new_am/2+1,new_am/2+1,new_am/2+1);
  fprintf(libderiv_header,"void init_libderiv_base();\n\n");
  fprintf(libderiv_header,"int  init_libderiv(Libderiv_t *, int max_am, int max_num_prim_quartets, int max_cart_class_size);\n");
  fprintf(libderiv_header,"void free_libderiv(Libderiv_t *);\n\n");
  fprintf(libderiv_header,"int  libderiv_storage_required(int max_am, int max_num_prim_quartets, int max_cart_class_size);\n");
  fprintf(libderiv_header,"#ifdef __cplusplus\n");
  fprintf(libderiv_header,"}\n");
  fprintf(libderiv_header,"#endif\n\n");
  fprintf(libderiv_header,"#endif\n");
  fclose(libderiv_header);
  ip_done();
  fclose(infile);
  fclose(outfile);
  exit(0);
}


void punt(char* str)
{
  printf(str);
  exit(1);
}


