/*
** fndcor(): C translation of the Fortran version, to remove the need to
**    link the library alloc, which also requires the linking of libparse,
**    etc, etc...
**
**
** This routine looks for the MEMORY keyword using the new-style input and
** libipv1.  
**
** Arguments: 
**    maxcrb  = int ptr to hold size of maxcore in bytes
**    infile  = file pointer to input file (eg input.dat)
**    outfile = file pointer to output file (eg output.dat)
**
** David Sherrill, February 1994
**
*/

#include <stdio.h>
#include <ip_libv1.h>

#define DEF_MAXCRR 100000    /* default maxcor in doubles: used only if
                              * it can't be read in */


void fndcor(maxcrb, infile, outfile)
      int *maxcrb ;
      FILE *infile, *outfile ;
{
char type[20] ;
char *s ;
int count ;
int maxcrr ;           /* maxcor in real words */
float size ;
int errcod ;
void fndcor_abort() ;

/*
   ffile(&infile,"input.dat",2);      * open for reading * 
   ffile(&outfile,"output.dat",1);    * open for append  * 
 
   ip_set_uppercase(1);
   ip_initialize(infile,outfile);
   ip_cwk_add(":DEFAULT");
   strcpy(progname, ":") ;
   strcat(progname, gprgid()) ;
   ip_cwk_add(progname) ;
*/

   maxcrr = DEF_MAXCRR ;  /* set maxcor to default first */

   errcod = ip_count("MEMORY", &count, 0) ;
   if (errcod == IPE_NOT_AN_ARRAY) {
      errcod = ip_data("MEMORY", "%d", &maxcrr, 0) ;
      if (errcod != IPE_OK) fndcor_abort(infile, outfile) ;
      }

   else if (errcod != IPE_OK) fndcor_abort(infile, outfile) ;

   else if (count == 1) {
      errcod = ip_data("MEMORY", "%d", &maxcrr, 1, 0) ;
      if (errcod != IPE_OK) fndcor_abort(infile, outfile) ;
      }

   else if (count == 2) {
      errcod = ip_data("MEMORY", "%f", &size, 1, 0) ;
      if (errcod != IPE_OK) fndcor_abort(infile, outfile) ;
      errcod = ip_data("MEMORY", "%s", type, 1, 1) ;
      if (errcod != IPE_OK) fndcor_abort(infile, outfile) ;
      /* convert string to uppercase */
      for (s=type; *s!='\0'; s++) {
         if (*s>='a' && *s <='z') *s = *s + 'A' - 'a';
         }
      if ((strcmp(type, "R")==0) || (strcmp(type, "REAL")==0))
         maxcrr = (int) size ;
      else if ((strcmp(type, "I")==0) || (strcmp(type, "INTEGER")==0)) 
         maxcrr = (int) (size * sizeof(int) / sizeof(double)) ;
      else if ((strcmp(type, "B")==0) || (strcmp(type, "BYTES")==0)) 
         maxcrr = (int) (size / sizeof(double)) ;
      else if ((strcmp(type, "KB")==0) || (strcmp(type, "KBYTES")==0)) 
         maxcrr = (int) (1000.0 * size / sizeof(double)) ;
      else if ((strcmp(type, "MB")==0) || (strcmp(type, "MBYTES")==0))
         maxcrr = (int) (1000000.0 * size / sizeof(double)) ;
      else {
         fprintf(outfile, "bad data type, specify one of: \n") ;
         fprintf(outfile, "REAL, INTEGER, BYTES, KBYTES, or MBYTES\n") ;
         fndcor_abort(infile, outfile) ;
         }
      }
  
   *maxcrb = maxcrr * sizeof(double) ;

/*   ip_done() ;
   fclose(infile) ;
   fclose(outfile) ;
*/
}


void fndcor_abort(infile, outfile)
      FILE *infile, *outfile ;
{
   fprintf(stderr, "Error: can't read MEMORY keyword!\n") ;
   fprintf(outfile, "Error: can't read MEMORY keyword!\n") ;
   ip_done() ;
   fclose(infile) ;
   fclose(outfile) ;
   exit(0) ;
}
 

