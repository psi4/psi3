/*
**
**  UGEOM: Read a geometry from geom.dat and write it to file30
**
*/

#include <stdio.h>
#include <stdlib.h>
#include <ip_libv1.h>
#include <libciomr.h>
#include <file30.h>

/* Function prototypes */
void init_io(void);
void exit_io(void);

FILE *infile, *outfile;

int main(int argc, char *argv[])
{
  int i,j;
  int natom, num, junk;
  double **geom, *zvals;
  char *geom_string;
  FILE *geometry;

  init_io();

  /* Open file30 and grab some info */
  file30_init();
  natom = file30_rd_natom();
  geom = file30_rd_geom();
  zvals = file30_rd_zvals();
  num = file30_rd_disp();
  file30_close();

  fprintf(outfile, "Old Geometry in File30:\n");
  for(i=0; i < natom; i++) {
      fprintf(outfile, "\n %1.0f ", zvals[i]);
      for(j=0; j < 3; j++)
          fprintf(outfile, "%20.10f  ", geom[i][j]);
    }
  fprintf(outfile, "\n\n");

  /* Check command line for desired geometry -- NB: C ordering */
  if(argc > 1) {
      if(strcmp(argv[1], "-n") == 0) {
         num = atoi(argv[2]);  
        }
    }
  geom_string = (char *) malloc(8 + (num > 9 ? 2 : 1));
  sprintf(geom_string, "%s%d", "GEOMETRY", num);

  /* Append geom.dat to the IP tree */
  ffile(&geometry, "geom.dat", 2);
  ip_append(geometry, outfile);

  /* Grab the desired geometry from geom.dat */
  ip_cwk_add(":OPTKING");
  if(!ip_exist(geom_string,0)) { 
     printf("No such geometry in geom.dat!\n");
     exit(2);
    }

  ip_count(geom_string, &junk, 0);
  if(junk != natom) {
     printf("Inconsistent atom numbers between geom.dat and file30!\n");
     exit(2);
    }
  for(i=0; i < natom; i++) {
      for(j=0; j < 3; j++)
          ip_data(geom_string, "%lf", &(geom[i][j]), 2, i, j);
    }

  /* Write the geometry to file30 and to output */
  file30_init();
  file30_wt_geom(geom);
  file30_wt_disp(++num);  /* Here we assume all is well! */
  file30_close();
  
  fprintf(outfile, "New Geometry from geom.dat:\n");
  for(i=0; i < natom; i++) {
      fprintf(outfile, "\n %1.0f ", zvals[i]);
      for(j=0; j < 3; j++)
          fprintf(outfile, "%20.10f  ", geom[i][j]);
    }
  fprintf(outfile, "\n\n");

  fprintf(outfile, "Next displacement: %d\n", num);

  /* Cleanup */
  fclose(geometry);
  free(geom_string);
  free(zvals);
  free_block(geom);

  exit_io();
  exit(0);
}

void init_io(void)
{
  int i;
  char *gprgid();
  char *progid;

  progid = (char *) malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s",gprgid());

  ffile(&infile,"input.dat",2);
  ffile(&outfile,"output.dat",1);
  tstart(outfile);
  ip_set_uppercase(1);
  ip_initialize(infile,outfile);
  ip_cwk_add(":DEFAULT");
  ip_cwk_add(progid);

  free(progid);
}

void exit_io(void)
{
  int i;
 
  ip_done();
  tstop(outfile);
  fclose(infile);
  fclose(outfile);
}

char *gprgid()
{
   char *prgid = "UGEOM";

   return(prgid);
}
