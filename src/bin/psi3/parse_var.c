/*
 * parse_var.c
 *
 * Here we try to parse the exec list from psi.dat or the user's exec
 */

#include <stdio.h>
#include <stdlib.h>
#include <libipv1/ip_lib.h>

#define MAX_VAR_SIZE 80

extern FILE *outfile;
extern void psi3_abort();

/*
 * This function will parse variable "name" and return an array
 * of variables that result.  The number of variables parsed is
 * returned to nvars, and mxvars gives the maximum number of vars.
 * This function needs to work in a recursive mode.
 */
char **parse_var(int *nvars, int mxvars, char *name)
{
  char **vars, **sub_vars;
  int num_vars=0, sub_num_vars, count, i, idx;
  int errcod;

  /* allocate enough space for an array of chars*/
  vars = (char **) malloc(mxvars*sizeof(char *)); 

  if (!ip_exist(name,0)) {
    fprintf(outfile,"Error (parse_var): couldn't locate %s.\n", name);
    psi3_abort();
  }

  errcod = ip_count(name,&count,0);
  if (errcod != IPE_OK && errcod != IPE_NOT_AN_ARRAY) {
    fprintf(outfile,"Error: trouble reading %s\n", name);
    psi3_abort();
  }
  if (errcod == IPE_NOT_AN_ARRAY) {
    /* vars[0] = (char *) malloc(MAX_VAR_SIZE*sizeof(char)); */
    /* ip_data(name,"%s",vars[0],0); */
    ip_string(name,&(vars[0]),0);
    if (vars[0][0] != '$') {
      *nvars = 1;
      return(vars);
    }
    else {
      sub_vars = parse_var(&num_vars, mxvars, &(vars[0][1]));
      if (num_vars > mxvars) {
        fprintf(outfile,"Error: num_vars = %d > mxvars\n", num_vars);
        psi3_abort();
      }
      free(vars[0]);
      for (i=0; i<num_vars; i++) {
        vars[i] = sub_vars[i];
      }
      free(sub_vars);
      *nvars = num_vars;
      return(vars);
    }
  } /* end the case of not an array */

  /* if we get this far, we must have an array */
  for (idx=0; idx<count; idx++) {
    /* vars[num_vars] = (char *) malloc(MAX_VAR_SIZE*sizeof(char)); */
    /* errcod = ip_data(name,"%s",vars[num_vars],1,idx); */
    errcod = ip_string(name,&(vars[num_vars]),1,idx);
    if (errcod != IPE_OK) {
      fprintf(outfile,"Error: can't read %s element %d\n", name, idx);
      psi3_abort();
    }
    if (vars[num_vars][0] != '$') num_vars++;
    else { 
      sub_vars = parse_var(&sub_num_vars, mxvars, &(vars[num_vars][1]));
      if (sub_num_vars+num_vars > mxvars) {
        fprintf(outfile,"Error: num_vars = %d > mxvars\n", 
                sub_num_vars+num_vars);
        psi3_abort();
      }
      free(vars[num_vars]);
      for (i=0; i<sub_num_vars; i++) {
        vars[num_vars+i] = sub_vars[i];
      }
      num_vars += sub_num_vars;
      free(sub_vars);
    }
  }  

  *nvars = num_vars;
  return(vars);
}



