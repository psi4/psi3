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

FILE *outfile, *vrr_header, *hrr_header, *libr12_header, *init_code;

void punt();
void emit_vrr_r_build(int,int,int);
void emit_vrr_t_build(int,int,int);
int emit_grt_order(int,int);
int emit_gr_order(int,int);
void emit_hrr_t_build(int,int);

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
  outfile = fopen("./output.dat", "w");
  hrr_header = fopen("./r12_hrr_header.h","w");
  vrr_header = fopen("./r12_vrr_header.h","w");
  libr12_header = fopen("./libr12.h","w");
  init_code = fopen("./init_libr12.c","w");

  /*---------------------------------
    Getting the new_am from libint.h
   ---------------------------------*/
  new_am = (MAX_AM - 1 - DERIV_LVL)*2;

  /* Setting up init_libr12.c, header.h */
  fprintf(init_code,"#include \"libr12.h\"\n");
  fprintf(init_code,"#include \"r12_hrr_header.h\"\n\n");
  fprintf(init_code,"/* This function initializes a matrix of pointers to routines */\n");
  fprintf(init_code,"/* for computing integral classes up to (%c%c|%c%c) */\n\n",
	  am_letter[new_am],am_letter[new_am],am_letter[new_am],am_letter[new_am]);
  fprintf(init_code,"void (*build_gr_abcd[%d][%d][%d][%d])();\n",new_am/2+1,new_am/2+1,new_am/2+1,new_am/2+1);
  fprintf(init_code,"void (*build_grt_abcd[%d][%d][%d][%d])();\n",new_am/2+1,new_am/2+1,new_am/2+1,new_am/2+1);
  fprintf(init_code,"prim_data *Data;\n");
  fprintf(init_code,"double *t1vrr_classes[%d][%d];\n",1+new_am,1+new_am);
  fprintf(init_code,"double *t2vrr_classes[%d][%d];\n",1+new_am,1+new_am);
  fprintf(init_code,"double *rvrr_classes[%d][%d];\n",1+new_am,1+new_am);
  fprintf(init_code,"double *gvrr_classes[%d][%d];\n",1+new_am,1+new_am);
  fprintf(init_code,"double *r12vrr_stack;\n\n");
  fprintf(init_code,"void init_libr12()\n{\n");

/*  emit_gr_order(0,new_am); */
  stack_size = emit_grt_order(0,new_am);
  emit_hrr_t_build(0,new_am);
  emit_vrr_r_build(0,new_am,max_class_size);
  emit_vrr_t1_build(0,new_am,max_class_size);
  emit_vrr_t2_build(0,new_am,max_class_size);
  
  fprintf(init_code,"}\n");
  fclose(init_code);
  fclose(hrr_header);
  fclose(vrr_header);
  
    /* Setting up libr12.h */
  fprintf(libr12_header,"/* Maximum angular momentum of functions in a basis set plus 1 */\n");
  fprintf(libr12_header,"#ifdef MAX_AM\n");
  fprintf(libr12_header," #undef MAX_AM\n");
  fprintf(libr12_header,"#endif\n");
  fprintf(libr12_header,"#define MAX_AM %d\n",1+new_am/2);
  fprintf(libr12_header,"#ifdef STACK_SIZE\n");
  fprintf(libr12_header," #undef STACK_SIZE\n");
  fprintf(libr12_header,"#endif\n");
  fprintf(libr12_header,"#define STACK_SIZE %d\n\n",stack_size);
  fprintf(libr12_header,"extern void (*build_gr_abcd[%d][%d][%d][%d])();\n",new_am/2+1,new_am/2+1,new_am/2+1,new_am/2+1);
  fprintf(libr12_header,"extern void (*build_grt_abcd[%d][%d][%d][%d])();\n",new_am/2+1,new_am/2+1,new_am/2+1,new_am/2+1);
  fprintf(libr12_header,"void init_libr12();\n\n");
  fclose(libr12_header);
  fclose(outfile);
  exit(0);
}


void punt(char* str)
{
  printf(str);
  exit(1);
}


