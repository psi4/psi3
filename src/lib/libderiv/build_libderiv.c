#include <math.h>
#include <stdio.h>
#include <ip_libv1.h>
#include <libciomr.h>
#include <libint.h>
#include "build_libderiv.h"

FILE *infile, *outfile, *vrr_header, *hrr_header, *dhrr_header,
     *deriv_header, *libderiv_header, *init_code;

void punt();
int emit_vrr_build(int,int,int);
int emit_vrr_order(int,int);
int emit_hrr_build(int,int);

int main()
{
  int i,j,k,l,f;
  int j_min, j_max, k_min, k_max, l_min, l_max;
  int errcod;
  int new_am, old_am;
  int class_size;
  int num_subfunctions;
  int max_class_size = 785;
  int stack_size;
  const int io[] = {0,1,3,6,10,15,21,28,36,45,55,66,78,91,105,120,136,153};
  const char am_letter[] = "0pdfghiklmnoqrtuvwxyz";


  /*-------------------------------
    Initialize files and libraries
   -------------------------------*/
  infile = fopen("./input.dat", "r");
  outfile = fopen("./output.dat", "w");
  dhrr_header = fopen("./dhrr_header.h","w");
  deriv_header = fopen("./deriv_header.h","w");
  libderiv_header = fopen("./libderiv.h","w");
  init_code = fopen("./init_libderiv.c","w");

  ip_set_uppercase(1);
  ip_initialize(infile,outfile);
  ip_cwk_add(":DEFAULT");
  ip_cwk_add(":LIBDERIV");
  
  /*----------------------------------------
    Getting the new_am from user and making
    sure it is consistent with libint.h
   ----------------------------------------*/
  errcod = ip_data("NEW_AM","%d",&new_am,0);
  if (errcod != IPE_OK)
    new_am = DEFAULT_NEW_AM;
  if (new_am <= 0)
    punt("  NEW_AM must be positive.");
  if (new_am > (LIBINT_MAX_AM - 1 - DERIV_LVL)*2)
    punt("  Maximum NEW_AM is greater installed libint.a allows.\n  Recompile libint.a with greater NEW_AM.");

  /* Setting up init_libderiv.c, header.h */
  fprintf(init_code,"#include <stdlib.h>\n");
  fprintf(init_code,"#include <string.h>\n");
  fprintf(init_code,"#include <libint.h>\n");
  fprintf(init_code,"#include \"libderiv.h\"\n");
  fprintf(init_code,"#include \"dhrr_header.h\"\n\n");
  fprintf(init_code,"/* This function initializes a matrix of pointers to routines */\n");
  fprintf(init_code,"/* for computing derivatives of ERI classes up to (%cs|%cs) */\n\n",am_letter[new_am],am_letter[new_am]);
  fprintf(init_code,"void (*build_deriv1_eri[%d][%d][%d][%d])(Libderiv_t *, int);\n",new_am/2+1,new_am/2+1,new_am/2+1,new_am/2+1);
  fprintf(init_code,"void init_libderiv_base()\n{\n");


  stack_size = emit_order(0,new_am);
  emit_dhrr_build(0,new_am);
  emit_deriv_build(0,new_am);

  
  fprintf(init_code,"}\n");
  fprintf(init_code,"/* These functions initialize library objects */\n");
  fprintf(init_code,"/* Library objects operate independently of each other */\n");
  fprintf(init_code,"int init_libderiv(Libderiv_t *libderiv, int num_prim_quartets, int max_cart_class_size)\n{\n");
  fprintf(init_code,"  int memory = 0;\n\n");
  fprintf(init_code,"  libderiv->int_stack = (double *) malloc(STACK_SIZE*sizeof(double));\n");
  fprintf(init_code,"  memory += STACK_SIZE;\n");
  fprintf(init_code,"  libderiv->zero_stack = (double *) malloc(max_cart_class_size*sizeof(double));\n");
  fprintf(init_code,"  bzero((char *)libderiv->zero_stack,max_cart_class_size*sizeof(double));\n");
  fprintf(init_code,"  memory += max_cart_class_size;\n");
  fprintf(init_code,"  libderiv->PrimQuartet = (prim_data *) malloc(num_prim_quartets*sizeof(prim_data));\n");
  fprintf(init_code,"  memory += num_prim_quartets*sizeof(prim_data)/sizeof(double);\n");
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
  fprintf(init_code,"  return;\n}\n");
  fclose(init_code);
  fclose(dhrr_header);
  fclose(deriv_header);
  
    /* Setting up libderiv.h */
  fprintf(libderiv_header,"/* Maximum angular momentum of functions in a basis set plus 1 */\n");
  fprintf(libderiv_header,"#define LIBDERIV_MAX_AM %d\n",1+new_am/2);
  fprintf(libderiv_header,"#ifdef STACK_SIZE\n");
  fprintf(libderiv_header," #undef STACK_SIZE\n");
  fprintf(libderiv_header,"#endif\n");
  fprintf(libderiv_header,"#define STACK_SIZE %d\n",stack_size);
  fprintf(libderiv_header,"#ifdef DERIV_LVL\n");
  fprintf(libderiv_header," #undef DERIV_LVL\n");
  fprintf(libderiv_header,"#endif\n");
  fprintf(libderiv_header,"#define DERIV_LVL %d\n\n",DERIV_LVL);
  fprintf(libderiv_header,"typedef struct {\n");
  fprintf(libderiv_header,"  double *int_stack;\n"); 
  fprintf(libderiv_header,"  prim_data *PrimQuartet;\n");
  fprintf(libderiv_header,"  double *zero_stack;\n");
  fprintf(libderiv_header,"  double *ABCD[12];\n");
  fprintf(libderiv_header,"  double AB[3];\n");
  fprintf(libderiv_header,"  double CD[3];\n");
  fprintf(libderiv_header,"  double *deriv_classes[%d][%d][%d];\n",1+new_am,1+new_am,12);
  fprintf(libderiv_header,"  double *dvrr_classes[%d][%d];\n",1+new_am,1+new_am);
  fprintf(libderiv_header,"  double *dvrr_stack;\n");
  fprintf(libderiv_header,"  } Libderiv_t;\n\n");
  fprintf(libderiv_header,"extern void (*build_deriv1_eri[%d][%d][%d][%d])(Libderiv_t *, int);\n",
	  new_am/2+1,new_am/2+1,new_am/2+1,new_am/2+1);
  fprintf(libderiv_header,"void init_libderiv_base();\n\n");
  fprintf(libderiv_header,"int  init_libderiv(Libderiv_t *, int, int);\n");
  fprintf(libderiv_header,"void free_libderiv(Libderiv_t *);\n\n");
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


