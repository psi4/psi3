#include <stdio.h>
#include <math.h>
#define EXTERN
#include "globals.h"

double norm_C(dpdfile2 *CME, dpdfile2 *Cme,
    dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf)
{
  double norm = 0.0;

  norm += dpd_file2_dot_self(CME);
  norm += dpd_file2_dot_self(Cme);
  norm += dpd_buf4_dot_self(CMNEF);
  norm += dpd_buf4_dot_self(Cmnef);
  norm += dpd_buf4_dot_self(CMnEf);

  norm = sqrt(norm);
  return norm;
}

double norm_C_rhf(dpdfile2 *CME, dpdbuf4 *CMnEf, dpdbuf4 *CMnfE) {
  double norm = 0.0;
  norm = 2.0 * dpd_file2_dot_self(CME);
  norm += 2.0 * dpd_buf4_dot_self(CMnEf);
  norm -= dpd_buf4_dot(CMnEf, CMnfE);
  norm = sqrt(norm);
  return norm;
}

double norm_C1(dpdfile2 *CME, dpdfile2 *Cme)
{
  double norm = 0.0;

  norm += dpd_file2_dot_self(CME);
  norm += dpd_file2_dot_self(Cme);
  norm = sqrt(norm);

  return norm;
}

double norm_C1_rhf(dpdfile2 *CME)
{
  double norm = 0.0;

  norm = 2*dpd_file2_dot_self(CME);
  norm = sqrt(norm);

  return norm;
}

void scm_C(dpdfile2 *CME, dpdfile2 *Cme, dpdbuf4 *CMNEF,
    dpdbuf4 *Cmnef, dpdbuf4 *CMnEf, double a)
{
  dpd_file2_scm(CME,a);
  dpd_file2_scm(Cme,a);
  dpd_buf4_scm(CMNEF,a);
  dpd_buf4_scm(Cmnef,a);
  dpd_buf4_scm(CMnEf,a);
  return;
}

void scm_C2(dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf, double a)
{
  dpd_buf4_scm(CMNEF,a);
  dpd_buf4_scm(Cmnef,a);
  dpd_buf4_scm(CMnEf,a);
  return;
}

void scm_C1(dpdfile2 *CME, dpdfile2 *Cme, double a)
{
  dpd_file2_scm(CME,a);
  dpd_file2_scm(Cme,a);
  return;
}

