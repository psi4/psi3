/*!
  \file rd_lagr.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*!
** file30_rd_lagr():	Reads in a lagrangian matrix in MO basis
**
**   STB(10/29/99)-This is now a wrapper function for rd_alpha_lagr() since in cscf
**   if there is no beta then alpha acts as the restricted solution.
**
**   returns: 
**	double **lagr	a matrix nmo by nmo.
*/


double **file30_rd_lagr(void)
{
    return file30_rd_alpha_lagr();
}
