/*---------------------------------------------------------------------------------------------
  This code generates Level 0 and 1 routines for computing non-standard two-electron integrals
  necessary in linear R12-methods. Evaluation is performed in Cartesian environment
  (as opposed to Hermite).

  The following functions are performed here:
  1) emit_grt_order builds Level 1 routines which compute ERIs(g), r12-integrals(r), and
     integrals over [r12,T1] and [r12,T2] (t1 and t2).
  2) emit_gr_order build Level 1 routines which compute ERIs(g) and r12-integrals(r) only.
  3) emit_vrr_r_build builds Level 0 routines which compute "VRR" r-integrals (a0||b0)
  4) emit_vrr_t_build builds Level 0 routines which compute "VRR" t1- and t2- integrals
     (a0|[r12,T1]|b0) and (a0|[r12,T2]|b0).
  5) emit_hrr_t_build builds Level 0 routines which compute "HRR" t1- and t2-integrals
     (ab|[r12,T1]|cd) and (ab|[r12,T2]|cd).

  Regular VRR and HRR Level 0 routines reside in LIBINT

 ---------------------------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <ip_libv1.h>
#include <libciomr.h>
#include <libint.h>
#include "build_libr12.h"

FILE *infile, *outfile, *vrr_header, *hrr_header, *libr12_header, *init_code;

void punt();
void emit_vrr_r_build(int,int,int);
void emit_vrr_t_build(int,int,int);
int emit_grt_order(int,int,int);
int emit_gr_order(int,int,int);
void emit_hrr_t_build(int,int);

int main()
{
  int i,j,k,l,f;
  int j_min, j_max, k_min, k_max, l_min, l_max;
  int errcod;
  int new_am, opt_am;
  int class_size;
  int num_subfunctions;
  int max_class_size = DEFAULT_MAX_CLASS_SIZE;
  int stack_size;
  const int io[] = {0,1,3,6,10,15,21,28,36,45,55,66,78,91,105,120,136,153};
  const char am_letter[] = "0pdfghiklmnoqrtuvwxyz";


  /*-------------------------------
    Initialize files and libraries
   -------------------------------*/
  infile = fopen("./input.dat", "r");
  outfile = fopen("./output.dat", "w");
  hrr_header = fopen("./r12_hrr_header.h","w");
  vrr_header = fopen("./r12_vrr_header.h","w");
  libr12_header = fopen("./libr12.h","w");
  init_code = fopen("./init_libr12.c","w");

  ip_set_uppercase(1);
  ip_initialize(infile,outfile);
  ip_cwk_add(":DEFAULT");
  ip_cwk_add(":LIBR12");

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

  errcod = ip_data("OPT_AM","%d",&opt_am,0);
  if (errcod != IPE_OK || opt_am < 2)
    opt_am = DEFAULT_OPT_AM;
  if (opt_am > new_am) opt_am = new_am;

  errcod = ip_data("MAX_CLASS_SIZE","%d",&max_class_size,0);
  if (max_class_size < 10)
    punt("  MAX_CLASS_SIZE cannot be smaller than 10.");
  
  /* Setting up init_libr12.c, header.h */
  fprintf(init_code,"#include <stdlib.h>\n");
  fprintf(init_code,"#include <libint.h>\n");
  fprintf(init_code,"#include \"libr12.h\"\n");
  fprintf(init_code,"#include \"r12_hrr_header.h\"\n\n");
  fprintf(init_code,"/* This function initializes a matrix of pointers to routines */\n");
  fprintf(init_code,"/* for computing integral classes up to (%c%c|%c%c) */\n\n",
	  am_letter[new_am],am_letter[new_am],am_letter[new_am],am_letter[new_am]);
  fprintf(init_code,"void (*build_r12_gr[%d][%d][%d][%d])(Libr12_t *, int);\n",new_am/2+1,new_am/2+1,new_am/2+1,new_am/2+1);
  fprintf(init_code,"void (*build_r12_grt[%d][%d][%d][%d])(Libr12_t *, int);\n",new_am/2+1,new_am/2+1,new_am/2+1,new_am/2+1);
  fprintf(init_code,"void init_libr12_base()\n{\n");

  /* Declare generic build routines */
  fprintf(vrr_header,"REALTYPE *r_vrr_build_xxxx(int am[2], prim_data *, REALTYPE *, const REALTYPE *,");
  fprintf(vrr_header," const REALTYPE *, REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *);\n");
  fprintf(vrr_header,"REALTYPE *t1_vrr_build_xxxx(int am[2], prim_data *, contr_data *, REALTYPE *,");
  fprintf(vrr_header," REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *);\n");
  fprintf(vrr_header,"REALTYPE *t2_vrr_build_xxxx(int am[2], prim_data *, contr_data *, REALTYPE *,");
  fprintf(vrr_header," REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *);\n");
  
/*  emit_gr_order(0,new_am); */
  stack_size = emit_grt_order(0,new_am,opt_am);
  emit_hrr_t_build(new_am, max_class_size);
  /*--- VRR build routines are optimized for classes up to opt_am/2 ---*/
  emit_vrr_r_build(0,opt_am,max_class_size);
  emit_vrr_t1_build(0,opt_am,max_class_size);
  emit_vrr_t2_build(0,opt_am,max_class_size);

  fprintf(init_code,"}\n\n");
  fprintf(init_code,"/* These functions initialize library objects */\n");
  fprintf(init_code,"/* Library objects operate independently of each other */\n");
  fprintf(init_code,"int init_libr12(Libr12_t *libr12, int num_prim_quartets)\n{\n");
  fprintf(init_code,"  int memory = 0;\n\n");
  fprintf(init_code,"  libr12->int_stack = (REALTYPE *) malloc(STACK_SIZE*sizeof(REALTYPE));\n");
  fprintf(init_code,"  memory += STACK_SIZE;\n");
  fprintf(init_code,"  libr12->PrimQuartet = (prim_data *) malloc(num_prim_quartets*sizeof(prim_data));\n");
  fprintf(init_code,"  memory += num_prim_quartets*sizeof(prim_data)/sizeof(REALTYPE);\n");
  fprintf(init_code,"  return memory;\n}\n\n");
  fprintf(init_code,"void free_libr12(Libr12_t *libr12)\n{\n");
  fprintf(init_code,"  if (libr12->int_stack != NULL) {\n");
  fprintf(init_code,"    free(libr12->int_stack);\n");
  fprintf(init_code,"    libr12->int_stack = NULL;\n");
  fprintf(init_code,"  }\n");
  fprintf(init_code,"  if (libr12->PrimQuartet != NULL) {\n");
  fprintf(init_code,"    free(libr12->PrimQuartet);\n");
  fprintf(init_code,"    libr12->PrimQuartet = NULL;\n");
  fprintf(init_code,"  }\n\n");
  fprintf(init_code,"  return;\n}\n");
  fclose(init_code);
  fclose(hrr_header);
  fclose(vrr_header);
  
    /* Setting up libr12.h */
  fprintf(libr12_header,"/* Maximum angular momentum of functions in a basis set plus 1 */\n");
  fprintf(libr12_header,"#define LIBR12_MAX_AM %d\n",1+new_am/2);
  fprintf(libr12_header,"#define LIBR12_OPT_AM %d\n",1+opt_am/2);
  fprintf(libr12_header,"#ifdef STACK_SIZE\n");
  fprintf(libr12_header," #undef STACK_SIZE\n");
  fprintf(libr12_header,"#endif\n");
  fprintf(libr12_header,"#define STACK_SIZE %d\n",stack_size);
  fprintf(libr12_header,"#define NUM_TE_TYPES 4\n\n");
  fprintf(libr12_header,"typedef struct {\n");
  fprintf(libr12_header,"  REALTYPE AB[3];\n");
  fprintf(libr12_header,"  REALTYPE CD[3];\n");
  fprintf(libr12_header,"  REALTYPE AC[3];\n");
  fprintf(libr12_header,"  REALTYPE ABdotAC, CDdotCA;\n");
  fprintf(libr12_header,"  } contr_data;\n\n");
  fprintf(libr12_header,"typedef struct {\n");
  fprintf(libr12_header,"  REALTYPE *int_stack;\n"); 
  fprintf(libr12_header,"  prim_data *PrimQuartet;\n");
  fprintf(libr12_header,"  contr_data ShellQuartet;\n");
  fprintf(libr12_header,"  REALTYPE *te_ptr[NUM_TE_TYPES];\n");
  fprintf(libr12_header,"  REALTYPE *t1vrr_classes[%d][%d];\n",1+new_am,1+new_am);
  fprintf(libr12_header,"  REALTYPE *t2vrr_classes[%d][%d];\n",1+new_am,1+new_am);
  fprintf(libr12_header,"  REALTYPE *rvrr_classes[%d][%d];\n",1+new_am,1+new_am);
  fprintf(libr12_header,"  REALTYPE *gvrr_classes[%d][%d];\n",2+new_am,2+new_am);
  fprintf(libr12_header,"  REALTYPE *r12vrr_stack;\n\n");
  fprintf(libr12_header,"  } Libr12_t;\n\n");
  fprintf(libr12_header,"extern void (*build_r12_gr[%d][%d][%d][%d])(Libr12_t *, int);\n",new_am/2+1,new_am/2+1,new_am/2+1,new_am/2+1);
  fprintf(libr12_header,"extern void (*build_r12_grt[%d][%d][%d][%d])(Libr12_t *, int);\n",new_am/2+1,new_am/2+1,new_am/2+1,new_am/2+1);
  fprintf(libr12_header,"void init_libr12_base();\n\n");
  fprintf(libr12_header,"int  init_libr12(Libr12_t *, int);\n");
  fprintf(libr12_header,"void free_libr12(Libr12_t *);\n\n");
  fclose(libr12_header);
  fclose(outfile);
  exit(0);
}


void punt(char* str)
{
  printf(str);
  exit(1);
}


