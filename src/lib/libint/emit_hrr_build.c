#include <math.h>
#include <stdio.h>
#include <libciomr.h>

#include "build_libint.h"

extern FILE *infile, *outfile, *hrr_header;

extern void punt(char *);
static int hash(int a[2][3], int b[2]);

int emit_hrr_build(int old_am, int new_am)
{
  FILE *code;
  int p,q,r,s;
  int ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz;
  int t0, t1, t2, t3, t4;
  int i,j,nj,i_i0,i_i1;
  int k,l,nl,k_i0,k_i1;
  int i0_step,i1_step;
  int a, b;
  int flag;
  int am_in[2];
  int am[2][3];
  int xyz;
  int class_size;
  int la, lb;
  int ld, lc, ld_max;
  int curr_count,curr_subfunction;
  static int io[] = {1,3,6,10,15,21,28,36,45,55,66,78,91,105,120,136,153};
  static const char am_letter[] = "0pdfghiklmnoqrtuvwxyz";
  char code_name[20];
  char function_name[18];
  

  for(lc=0;lc<=new_am;lc++) {
    ld_max = (lc+1)/2;
    for(ld=1;ld<=ld_max;ld++) {

      /*-----------------------
	HRR on centers C and D
       -----------------------*/

      am_in[0] = lc-ld;
      am_in[1] = ld;

      sprintf(function_name,"hrr3_build_%c%c",am_letter[am_in[0]],am_letter[am_in[1]]);
      sprintf(code_name,"hrr3_build_%c%c.c",am_letter[am_in[0]],am_letter[am_in[1]]);
      code = fopen(code_name,"w");

      /* include the function into the hrr_header.h */
      fprintf(hrr_header,"void %s(const double *, double *, const double *, const double *, int);\n",
	      function_name);

      fprintf(code,"  /* This machine-generated function computes a quartet of |%c%c) integrals */\n\n",
	      am_letter[am_in[0]],am_letter[am_in[1]]);

      fprintf(code,"void %s(const double *CD, double *vp, const double *I0, const double *I1, int ab_num)\n{\n",function_name);
      fprintf(code,"  int ab;\n");
      fprintf(code,"  const double CD0 = CD[0];\n");
      fprintf(code,"  const double CD1 = CD[1];\n");
      fprintf(code,"  const double CD2 = CD[2];\n");

      nl = (am_in[1]*(am_in[1]+1))/2;
      i0_step = (am_in[0]+2)*(am_in[0]+3)*nl/2;
      i1_step = (am_in[0]+1)*(am_in[0]+2)*nl/2;
      fprintf(code,"  for(ab=0;ab<ab_num;ab++) {\n");

      for(p = 0; p <= am_in[0]; p++){
	am[0][0] = am_in[0] - p;
	for(q = 0; q <= p; q++){
	  am[0][1] = p - q;
	  am[0][2] = q;
	  
	  for(r = 0; r <= am_in[1]; r++){
	    am[1][0] = am_in[1] - r;
	    for(s = 0; s <= r; s++){
	      am[1][1] = r - s;
	      am[1][2] = s;

	      if (am[1][0]) /* build along x */
		xyz = 0;
	      else if (am[1][1]) /* build along y */
		xyz = 1;
	      else /*build along z */
		xyz = 2;

	      am[0][xyz] += 1;
	      am_in[0] += 1;
	      am[1][xyz] -= 1;
	      am_in[1] -= 1;
	      t0 = hash(am,am_in);
	      am[0][xyz] -= 1;
	      am_in[0] -= 1;
	      t1 = hash(am,am_in);
	      am[1][xyz] += 1;
	      am_in[1] += 1;
	      
	      fprintf(code, "    *(vp++) = I0[%d] + CD%d*I1[%d];\n",t0,xyz,t1);
	    }
	  }
	}
      }
      fprintf(code,"    I0 += %d;\n    I1 += %d;\n",i0_step,i1_step);
      fprintf(code,"  }\n}\n");
      fclose(code);
      printf("Done with %s\n",code_name);

      
      
      /*-----------------------
	HRR on centers A and B
       -----------------------*/

      la = lc-ld;  lb = ld;
      am_in[0] = la;
      am_in[1] = lb;

      /* include the function into the hrr_header.h */
      fprintf(hrr_header,"void %s(const double *, double *, const double *, const double *, int);\n",
	      function_name);
      
      sprintf(function_name,"hrr1_build_%c%c",am_letter[am_in[0]],am_letter[am_in[1]]);
      sprintf(code_name,"hrr1_build_%c%c.c",am_letter[am_in[0]],am_letter[am_in[1]]);
      code = fopen(code_name,"w");

      fprintf(code,"  /* This machine-generated function computes a quartet of (%c%c| integrals */\n\n",
	      am_letter[am_in[0]],am_letter[am_in[1]]);
      fprintf(code,"void %s(const double *AB, double *vp, const double *I0, const double *I1, int cd_num)\n{\n",function_name);
      fprintf(code,"  int cd;\n");
      fprintf(code,"  const double *i0, *i1;\n");
      fprintf(code,"  const double AB0 = AB[0];\n");
      fprintf(code,"  const double AB1 = AB[1];\n");
      fprintf(code,"  const double AB2 = AB[2];\n");

      nj = (lb*(lb+1))/2;

      for(p = 0; p <= am_in[0]; p++){
	am[0][0] = am_in[0] - p;
	for(q = 0; q <= p; q++){
	  am[0][1] = p - q;
	  am[0][2] = q;
	  
	  for(r = 0; r <= am_in[1]; r++){
	    am[1][0] = am_in[1] - r;
	    for(s = 0; s <= r; s++){
	      am[1][1] = r - s;
	      am[1][2] = s;

	      if (am[1][0]) /* build along x */
		xyz = 0;
	      else if (am[1][1]) /* build along y */
		xyz = 1;
	      else /* build along z */
		xyz = 2;

	      am[0][xyz] += 1;
	      am_in[0] += 1;
	      am[1][xyz] -= 1;
	      am_in[1] -= 1;
	      t0 = hash(am,am_in);
	      am[0][xyz] -= 1;
	      am_in[0] -= 1;
	      t1 = hash(am,am_in);
	      am[1][xyz] += 1;
	      am_in[1] += 1;
	      
	      if (t0)
		fprintf(code,"  i0 = I0 + %d*cd_num;\n",t0);
	      else
		fprintf(code,"  i0 = I0;\n");
	      if (t1)
		fprintf(code,"  i1 = I1 + %d*cd_num;\n",t1);
	      else
		fprintf(code,"  i1 = I1;\n");

	      fprintf(code,"  for(cd=0;cd<cd_num;cd++)\n");
	      fprintf(code,"    *(vp++) = *(i0++) + AB%d*(*(i1++));\n",xyz);

	    }
	  }
	}
      }
      fprintf(code,"}\n");
      fclose(code);
      printf("Done with %s\n",code_name);
    }
  }
}


/*----------------------------------------------------------------------------------
  hash(a,b) returns a number of the (a[0] a[1]) type product within a doublet.
  a contains x y and z exponents of functions on centers A and B, and b contains
  their angular momenta
 ----------------------------------------------------------------------------------*/

int hash(a, b)
  int a[2][3];
  int b[2];
{
  int c[2] = {0,0};
  int i;
  static int io[] = {0,1,3,6,10,15,21,28,36,45,55,66,78,91,105,120,136,153};

  if(b[0]){
    i=b[0]-a[0][0];
    c[0]=i+io[i]-a[0][1];
    }
  if(b[1]){
    i=b[1]-a[1][0];
    c[1]=i+io[i]-a[1][1];
    }

  return c[0]*io[b[1]+1]+c[1];
}
