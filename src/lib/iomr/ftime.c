#include <stdio.h>
#include <libciomr.h>
void ftstart_()
   {
       FILE *outfile;
       ffile__(&outfile,"output.dat", 1);
       tstart__(outfile);
       fclose(outfile);
   }

void ftstop_()
   {
       FILE *outfile;
       ffile__(&outfile,"output.dat", 1);
       tstop__(outfile);
       fclose(outfile);
   }
void ftstart()
   {
       FILE *outfile;
       ffile__(&outfile,"output.dat", 1);
       tstart__(outfile);
       fclose(outfile);
   }

void ftstop()
   {
       FILE *outfile;
       ffile__(&outfile,"output.dat", 1);
       tstop__(outfile);
       fclose(outfile);
   }
