/*!
  \file fndcor.c
*/

#include <stdio.h>
#include <libipv1/ip_lib.h>

#define DEF_MAXCRR 100000    /* default maxcor in doubles: used only if
                              * it can't be read in */

static void fndcor_abort();

/*!
** fndcor(): C translation of the Fortran version, to remove the need to
**    link the library alloc, which also requires the linking of libparse,
**    etc, etc...
**
**
** This routine looks for the MEMORY keyword using the new-style input and
** libipv1.  
**
** Arguments: 
**    \param maxcrb  = long int ptr to hold size of maxcore in bytes
**    \param infile  = file pointer to input file (eg input.dat)
**    \param outfile = file pointer to output file (eg output.dat)
**
** David Sherrill, February 1994
** Revised to handle more than 2GB of memory by Ed Valeev, October 2000
**
*/

void fndcor(long int *maxcrb, FILE *infile, FILE *outfile)
{
  char type[20] ;
  char *s;
  int count ;
  long int maxcrr ;           /* maxcor in real words */
  char *maxcrr_str;           /* string representation of maxcrr */
  double size ;
  int errcod ;

   maxcrr = DEF_MAXCRR ;  /* set maxcor to default first */

   errcod = ip_count("MEMORY", &count, 0) ;
   if (errcod != IPE_OK) fndcor_abort(infile, outfile) ;
   else if (errcod == IPE_NOT_AN_ARRAY) { /* Scalar specification of MEMORY */
      errcod = ip_string("MEMORY", &maxcrr_str, 0);
      if (errcod != IPE_OK) fndcor_abort(infile, outfile) ;
      maxcrr = atol(maxcrr_str);
   }
   /* Array specification of MEMORY */
   else if (count == 1) {
      errcod = ip_string("MEMORY", &maxcrr_str, 0);
      if (errcod != IPE_OK) fndcor_abort(infile, outfile) ;
      maxcrr = atol(maxcrr_str);
   }
   else if (count == 2) {
      errcod = ip_data("MEMORY", "%lf", &size, 1, 0) ;
      if (errcod != IPE_OK) fndcor_abort(infile, outfile) ;
      errcod = ip_data("MEMORY", "%s", type, 1, 1) ;
      if (errcod != IPE_OK) fndcor_abort(infile, outfile) ;
      /* convert string to uppercase */
      for (s=type; *s!='\0'; s++) {
         if (*s>='a' && *s <='z') *s = *s + 'A' - 'a';
      }
      if ((strcmp(type, "R")==0) || (strcmp(type, "REAL")==0))
         maxcrr = (long int) size ;
      else if ((strcmp(type, "I")==0) || (strcmp(type, "INTEGER")==0)) 
         maxcrr = (long int) (size * sizeof(int) / sizeof(double)) ;
      else if ((strcmp(type, "B")==0) || (strcmp(type, "BYTES")==0)) 
         maxcrr = (long int) (size / sizeof(double)) ;
      else if ((strcmp(type, "KB")==0) || (strcmp(type, "KBYTES")==0)) 
         maxcrr = (long int) (1000.0 * size / sizeof(double)) ;
      else if ((strcmp(type, "MB")==0) || (strcmp(type, "MBYTES")==0))
         maxcrr = (long int) (1000000.0 * size / sizeof(double)) ;
      else if ((strcmp(type, "GB")==0) || (strcmp(type, "GBYTES")==0))
         maxcrr = (long int) (1000000000.0 * size / sizeof(double)) ;
      else {
         fprintf(outfile, "bad data type, specify one of: \n") ;
         fprintf(outfile, "REAL, INTEGER, BYTES, KBYTES, MBYTES, or GBYTES\n") ;
         fndcor_abort(infile, outfile) ;
         }
      }
  
   *maxcrb = maxcrr * sizeof(double) ;

   return;
}


static void fndcor_abort(infile, outfile)
      FILE *infile, *outfile ;
{
   fprintf(stderr, "Error: can't read MEMORY keyword!\n") ;
   fprintf(outfile, "Error: can't read MEMORY keyword!\n") ;
   ip_done() ;
   fclose(infile) ;
   fclose(outfile) ;
   exit(0) ;
}
 

