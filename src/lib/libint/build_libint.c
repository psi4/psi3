/*------------------------------------------------------------------------------------------------------
                                           BUILD_LIBINT
  		       Written by Dr. Justin T. Fermann and Edward F. Valeev
  This program generates files necessary for compiling the LIBINT library. LIBINT is a library of
  streamlined highly efficient routines for recursive computation of ERI integrals of the form (a0|b0).
  The library is used by the integral code CINTS. The following files must be present in the directory
  where the current program is intended to run: input.dat and Makefile.libint (all included with the source).
  input.dat contains set of input data necessary for building the library.
  LIBINT section of input.dat must contain the following keywords:
  NEW_AM - twice the desired maximum angular momentum (default is 8);
  OPT_AM - twice the angular momentum up to which VRR Level 0 routines are machine generated.
  A generic VRR Level 0 function is used past this value. OPT_AM=8 should be enough for almost anyone.
  MAX_CLASS_SIZE - maximum length of _build_a0b0 (if (a0|b0) class is longer than MAX_CLASS_SIZE, the routine is
  going to be split into several smaller ones. This is done to prevent compiler from exhausting system resources.
  Defaults to 785).

  *EXAMPLE*: if one wants LIBINT work for up to (gg|gg) integrals (the current PSI limit), NEW_AM has to be set
  to 8 (g corresponds to l=4, (gg|gg) class will require at most the (l0|l0) class). The intended angular
  momentum limit for PSI 3 is i (l=6), therefore up to (q0|q0) classes are required. NEW_AM must be set to 12.

  Accessing functions in LIBINT is very simple - the program has to call init_libint() just once before it
  starts computing integrals. After that all top_build_... functions (****BUT**** top_build_0000, which should
  never be neccessary, since (ss|ss) class is easily computed from the auxiliary function) will be arranged
  in a matrix of pointers. E.g., to call top_build_i0f0(...) one has to invoke top_build_a0b0[6][3](...).
 ------------------------------------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <ip_libv1.h>
#include <libciomr.h>
#include "build_libint.h"

FILE *infile, *outfile, *vrr_header, *hrr_header, *libint_header, *init_code;

void punt();
int emit_vrr_build(int,int,int);
int emit_vrr_order(int,int,int);
int emit_hrr_build(int,int);

int main()
{
  int i,j,k,l,f;
  int j_min, j_max, k_min, k_max, l_min, l_max;
  int errcod;
  int old_am = 0;
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
  vrr_header = fopen("./vrr_header.h","w");
  hrr_header = fopen("./hrr_header.h","w");
  libint_header = fopen("./libint.h","w");
  init_code = fopen("./init_libint.c","w");

  ip_set_uppercase(1);
  ip_initialize(infile,outfile);
  ip_cwk_add(":DEFAULT");
  ip_cwk_add(":LIBINT");

  /*------------------
    Parsing the input
   ------------------*/
  errcod = ip_data("NEW_AM","%d",&new_am,0);
  if (errcod != IPE_OK)
    new_am = DEFAULT_NEW_AM;
  if (new_am <= 0)
    punt("  NEW_AM must be positive.");
  if (new_am > MAX_AM)
    punt("  Maximum NEW_AM exceeded. Contact the program author.");

  errcod = ip_data("OPT_AM","%d",&opt_am,0);
  if (errcod != IPE_OK || opt_am < 2)
    opt_am = DEFAULT_OPT_AM;
  if (opt_am > new_am) opt_am = new_am;

  errcod = ip_data("MAX_CLASS_SIZE","%d",&max_class_size,0);
  if (max_class_size < 10)
    punt("  MAX_CLASS_SIZE cannot be smaller than 10.");

  /* Setting up init_libint.c, header.h */
  fprintf(init_code,"#include <stdio.h>\n");
  fprintf(init_code,"#include <stdlib.h>\n");
  fprintf(init_code,"#include \"libint.h\"\n");
  fprintf(init_code,"#include \"hrr_header.h\"\n\n");
  fprintf(init_code,"/* This function initializes a matrix of pointers to routines */\n");
  fprintf(init_code,"/* for computing ERI classes up to (%cs|%cs) - the base of the library */\n\n",
	  am_letter[new_am],am_letter[new_am]);
  fprintf(init_code,"double *(*build_eri[%d][%d][%d][%d])(Libint_t *, int);\n",new_am/2+1,new_am/2+1,new_am/2+1,new_am/2+1);
  fprintf(init_code,"void init_libint_base()\n{\n");

  /* Declare generic build routines */
  fprintf(vrr_header,"double *vrr_build_xxxx(int am[2], prim_data *, double *,const double *, const double *, const double *,
const double *, const double *);\n");

  stack_size = emit_order(old_am,new_am,opt_am);
  emit_vrr_build(old_am, opt_am, max_class_size);
  emit_hrr_build(old_am, new_am);
  
  fprintf(init_code,"}\n");
  fprintf(init_code,"/* These functions initialize library objects */\n");
  fprintf(init_code,"/* Library objects operate independently of each other */\n");
  fprintf(init_code,"int init_libint(Libint_t *libint, int num_prim_quartets)\n{\n");
  fprintf(init_code,"  int memory = 0;\n\n");
  fprintf(init_code,"  libint->int_stack = (double *) malloc(STACK_SIZE*sizeof(double));\n");
  fprintf(init_code,"  memory += STACK_SIZE;\n");
  fprintf(init_code,"  libint->PrimQuartet = (prim_data *) malloc(num_prim_quartets*sizeof(prim_data));\n");
  fprintf(init_code,"  memory += num_prim_quartets*sizeof(prim_data)/sizeof(double);\n");
  fprintf(init_code,"  return memory;\n}\n\n");
  fprintf(init_code,"void free_libint(Libint_t *libint)\n{\n");
  fprintf(init_code,"  if (libint->int_stack != NULL) {\n");
  fprintf(init_code,"    free(libint->int_stack);\n");
  fprintf(init_code,"    libint->int_stack = NULL;\n");
  fprintf(init_code,"  }\n");
  fprintf(init_code,"  if (libint->PrimQuartet != NULL) {\n");
  fprintf(init_code,"    free(libint->PrimQuartet);\n");
  fprintf(init_code,"    libint->PrimQuartet = NULL;\n");
  fprintf(init_code,"  }\n\n");
  fprintf(init_code,"  return;\n}\n");
  fclose(init_code);
  fclose(vrr_header);
  fclose(hrr_header);
  
    /* Setting up libint.h */
  fprintf(libint_header,"/* Maximum angular momentum of functions in a basis set plus 1 */\n");
  fprintf(libint_header,"#define LIBINT_MAX_AM %d\n",1+new_am/2);
  fprintf(libint_header,"#define LIBINT_OPT_AM %d\n",1+opt_am/2);
  fprintf(libint_header,"#define STACK_SIZE %d\n\n",stack_size);
  fprintf(libint_header,"typedef struct pdata{\n");
  fprintf(libint_header,"  double F[%d];\n",2*new_am+1);
  fprintf(libint_header,"  double U[6][3];\n");
  fprintf(libint_header,"  double twozeta_a;\n"); 
  fprintf(libint_header,"  double twozeta_b;\n"); 
  fprintf(libint_header,"  double twozeta_c;\n");
  fprintf(libint_header,"  double twozeta_d;\n");
  fprintf(libint_header,"  double oo2z;\n");
  fprintf(libint_header,"  double oo2n;\n");
  fprintf(libint_header,"  double oo2zn;\n");
  fprintf(libint_header,"  double poz;\n");
  fprintf(libint_header,"  double pon;\n");
  fprintf(libint_header,"  double oo2p;\n");
  fprintf(libint_header,"  double ss_r12_ss;\n");
  fprintf(libint_header,"  } prim_data;\n\n");
  fprintf(libint_header,"typedef struct {\n");
  fprintf(libint_header,"  double *int_stack;\n"); 
  fprintf(libint_header,"  prim_data *PrimQuartet;\n"); 
  fprintf(libint_header,"  double AB[3];\n");
  fprintf(libint_header,"  double CD[3];\n");
  fprintf(libint_header,"  double *vrr_classes[%d][%d];\n",1+new_am,1+new_am);
  fprintf(libint_header,"  double *vrr_stack;\n");
  fprintf(libint_header,"  } Libint_t;\n\n");
  fprintf(libint_header,"extern double *(*build_eri[%d][%d][%d][%d])(Libint_t *, int);\n",
	  new_am/2+1,new_am/2+1,new_am/2+1,new_am/2+1);
  fprintf(libint_header,"void init_libint_base();\n");
  fprintf(libint_header,"int  init_libint(Libint_t *, int);\n");
  fprintf(libint_header,"void free_libint(Libint_t *);\n\n");
  fclose(libint_header);
  ip_done();
  fclose(outfile);
  fclose(infile);
  exit(0);
}


void punt(char* str)
{
  printf(str);
  exit(1);
}


