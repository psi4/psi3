#include <math.h>
#include <stdio.h>
#include <ip_libv1.h>
#include <libciomr.h>
#include <libint.h>
#include "build_libderiv.h"

FILE *outfile, *vrr_header, *hrr_header, *dhrr_header, *deriv_header, *libderiv_header, *init_code;

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
  outfile = fopen("./output.dat", "w");
  dhrr_header = fopen("./dhrr_header.h","w");
  deriv_header = fopen("./deriv_header.h","w");
  libderiv_header = fopen("./libderiv.h","w");
  init_code = fopen("./init_libderiv.c","w");

  /*---------------------------------
    Getting the new_am from libint.h
   ---------------------------------*/
  new_am = (MAX_AM - 1 - DERIV_LVL)*2;

  /* Setting up init_libderiv.c, header.h */
  fprintf(init_code,"#include <libint.h>\n");
  fprintf(init_code,"#include \"libderiv.h\"\n");
  fprintf(init_code,"#include \"dhrr_header.h\"\n\n");
  fprintf(init_code,"/* This function initializes a matrix of pointers to routines */\n");
  fprintf(init_code,"/* for computing derivatives of ERI classes up to (%cs|%cs) */\n\n",am_letter[new_am],am_letter[new_am]);
  fprintf(init_code,"void (*dbuild_abcd[%d][%d][%d][%d])();\n",new_am/2+1,new_am/2+1,new_am/2+1,new_am/2+1);
  fprintf(init_code,"prim_data *Data;\n");
  fprintf(init_code,"double *deriv_classes[%d][%d][%d];\n",1+new_am,1+new_am,12);
  fprintf(init_code,"double *dvrr_classes[%d][%d];\n",1+new_am,1+new_am);
  fprintf(init_code,"double *dvrr_stack;\n\n");
  fprintf(init_code,"void init_libderiv()\n{\n");


  stack_size = emit_order(0,new_am);
  emit_dhrr_build(0,new_am);
  emit_deriv_build(0,new_am);

  
  fprintf(init_code,"}\n");
  fclose(init_code);
  fclose(dhrr_header);
  fclose(deriv_header);
  
    /* Setting up libderiv.h */
  fprintf(libderiv_header,"/* Maximum angular momentum of functions in a basis set plus 1 */\n");
/*  fprintf(libderiv_header,"#define MAX_AM %d\n",1+new_am/2);*/
  fprintf(libderiv_header,"#ifdef STACK_SIZE\n");
  fprintf(libderiv_header," #undef STACK_SIZE\n");
  fprintf(libderiv_header,"#endif\n");
  fprintf(libderiv_header,"#define STACK_SIZE %d\n",stack_size);
  fprintf(libderiv_header,"#ifdef DERIV_LVL\n");
  fprintf(libderiv_header," #undef DERIV_LVL\n");
  fprintf(libderiv_header,"#endif\n");
  fprintf(libderiv_header,"#define DERIV_LVL %d\n\n",DERIV_LVL);
/*  fprintf(libderiv_header,"typedef struct {\n");
  fprintf(libderiv_header,"  double F[%d];\n",4*MAX_AM-3);
  fprintf(libderiv_header,"  double U[6][3];\n");
  fprintf(libderiv_header,"  double twozeta_a;\n");
  fprintf(libderiv_header,"  double twozeta_b;\n");
  fprintf(libderiv_header,"  double twozeta_c;\n");
  fprintf(libderiv_header,"  double twozeta_d;\n");
  fprintf(libderiv_header,"  double oo2z;\n");
  fprintf(libderiv_header,"  double oo2n;\n");
  fprintf(libderiv_header,"  double oo2zn;\n");
  fprintf(libderiv_header,"  double poz;\n");
  fprintf(libderiv_header,"  double pon;\n");
  fprintf(libderiv_header,"  double oo2p;\n");
  fprintf(libderiv_header,"  double ss_r12_ss;\n");
  fprintf(libderiv_header,"  } prim_data;\n\n");*/
  fprintf(libderiv_header,"extern void (*dbuild_abcd[%d][%d][%d][%d])();\n",new_am/2+1,new_am/2+1,new_am/2+1,new_am/2+1);
  fprintf(libderiv_header,"void init_libderiv();\n\n");
  fclose(libderiv_header);
  fclose(outfile);
  exit(0);
}


void punt(char* str)
{
  printf(str);
  exit(1);
}


