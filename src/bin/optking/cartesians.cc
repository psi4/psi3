/******************************************************************************

       cartesians.cc

       this file contains member functions for cartesian not included
       in the declaration (cartesians.h)

******************************************************************************/ 


#include <math.h>

extern "C" {
  #include <stdio.h>
  #include <stdlib.h>
  #include <string.h>
  #include <libciomr/libciomr.h>
  #include <libipv1/ip_lib.h>
  #include <libfile30/file30.h>
  #include <physconst.h>
  #include <masses.h>
}

#define EXTERN
#include "opt.h"
#include "cartesians.h"





/*-----------------------------------------------------------------------------
       FORCES

       This function returns the forces in cartesian coordinates in aJ/Ang

-----------------------------------------------------------------------------*/ 

double *cartesians:: forces() {
  int i;
  double *f;
  f = init_array(3*num_atoms);

  for (i=0;i<3*num_atoms;++i)
    f[i] = -1.0 * grad[i] * _hartree2J * 1.0E18 / _bohr2angstroms;
  return f;
}





/*-----------------------------------------------------------------------------
       CARTESIANS this is the class constructor for cartesian
----------------------------------------------------------------------------*/
cartesians::cartesians() {
  int i, j, a, count, isotopes_given = 0, masses_given = 0;
  char label[MAX_LINELENGTH], line1[MAX_LINELENGTH];
  FILE *fp_11;
  double an,x,y,z,tval;
  double **geom, *zvals;

  if (   ((fp_11 = fopen("file11.dat","r")) == NULL)
      || (optinfo.numerical_dertype > 0) ) {
     /* Read geometry and energy data only from file30 */
    rewind(fp_input);
    ip_set_uppercase(1);
    ip_initialize(fp_input,outfile);
    ip_cwk_clear();
    ip_cwk_add(":DEFAULT");
    ip_cwk_add(":OPTKING");

    file30_init();
    geom = file30_rd_geom();
    num_atoms = file30_rd_natom();
    zvals = file30_rd_zvals();
/*** FIX ***/
    energy = file30_rd_escf();
    file30_close();
    ip_done();

     coord = new double[3*num_atoms];
     grad = new double[3*num_atoms];
     atomic_num = new double[num_atoms];
     mass = new double[3*num_atoms];
     for(i=0, count=-1; i < num_atoms; i++) {
         atomic_num[i] = zvals[i];
         for(j=0; j < 3; j++) {
             coord[++count] = geom[i][j];
           }
       }
     free(zvals);
     free_block(geom);
  }
  else {  /* Get what you need from file11 */

  if (fgets(label, MAX_LINELENGTH, fp_11) == NULL) {
     fprintf(outfile, "Trouble reading first line of file11.dat\n");
     exit(0);
  }

  /* now read the number of atoms and the energy */
  fgets(line1, MAX_LINELENGTH, fp_11);
  if (sscanf(line1, "%d %lf", &num_atoms, &energy) != 2) {
     fprintf(outfile,
       "Trouble reading natoms and energy from second line of file11.dat\n");
     exit(0);
  }

  coord = new double[3*num_atoms];
  grad = new double[3*num_atoms];
  mass = new double[3*num_atoms];
  atomic_num = new double[num_atoms];

  /* now rewind file11 */
  rewind(fp_11) ;

  /* read in one chunk at a time */
   while (fgets(label, MAX_LINELENGTH, fp_11) != NULL) {
      fgets(line1, MAX_LINELENGTH, fp_11) ;
      if (sscanf(line1, "%d %lf", &num_atoms, &energy) != 2) {
         fprintf(outfile,
         "(cartesians()): Trouble reading natoms and energy from file11.dat\n");
         exit(0);
      }
      count = -1;
      for (i=0; i<num_atoms ; i++) {
         if (fscanf(fp_11, "%lf %lf %lf %lf", &an, &x, &y, &z) != 4) {
            fprintf(outfile,
              "(cartesians()): Trouble reading coordinates from file11.dat\n");
            exit(0);
         }
         atomic_num[i]=an; coord[++count]=x; coord[++count]=y; coord[++count]=z;
      }
      count = -1;
      for (i=0; i<num_atoms ; i++) {
         if (fscanf(fp_11, "%lf %lf %lf", &x, &y, &z) != 3) {
            fprintf(outfile,
              "(cartesians()): Trouble reading gradients from file11.dat\n");
            exit(0);
         }
         grad[++count] = x; grad[++count] = y; grad[++count] = z;
      }
      fgets(line1, MAX_LINELENGTH, fp_11);
   } 
   fclose(fp_11) ;

   } /* else */

  /* read masses from input.dat or use default masses */
   count = -1;
   rewind(fp_input);
   ip_set_uppercase(1);
   ip_initialize(fp_input,outfile);
   ip_cwk_clear();
   ip_cwk_add(":DEFAULT");
   ip_cwk_add(":OPTKING");
   if (ip_exist("ISOTOPES",0)) {
      isotopes_given = 1;
      a = 0;
      ip_count("ISOTOPES", &a, 0);
      if (a != num_atoms) {
         fprintf(outfile,"ISOTOPES array has wrong dimension.\n");
         exit(2);
      }
      for (i=0;i<num_atoms;++i) {
         ip_data("ISOTOPES","%s", line1,1,i);
         for (j=0;j<138;j++) {
            if (strcmp(line1, mass_labels[j]) == 0) {
               mass[++count] = atomic_masses[j];
               mass[++count] = atomic_masses[j];
               mass[++count] = atomic_masses[j];
               break;
            }
         }
         if (j == 138) {
            fprintf(outfile,
              "Isotope label %s is unidentifiable.\n",mass_labels[j]);
            exit(2);
         }
      }
   }
   if (ip_exist("MASSES",0)) {
      masses_given = 1;
      if (isotopes_given)
         fprintf(outfile,"Ignoring ISOTOPES keyword, using given MASSES.\n");
      a = 0;
      ip_count("MASSES",&a,0);
      if (a != num_atoms) {
         fprintf(outfile,"MASSES array has wrong dimension\n");
         exit(2);
      }
      else {
         for(i=0;i<num_atoms;++i) {
            ip_data("MASSES","%lf",&tval,1,i);
            mass[++count] = tval;
            mass[++count] = tval;
            mass[++count] = tval;
         }
      }
   }
   if ((isotopes_given == 0) && (masses_given == 0)) {
      for(i=0;i<num_atoms;++i) {
         a = (int) get_atomic_num(i); // casting to an int for index
         mass[++count] = an2masses[a];
         mass[++count] = an2masses[a];
         mass[++count] = an2masses[a];
      }
   }
   ip_done();
   return;
}





/*-----------------------------------------------------------------------------

       PRINT

       flag = 0 print geometry to fp_out
       flag = 1 print Z, geom to fp_out
       flag = 2 print geom and grad with masses
       flag = 4 print geometry to geom.dat
       flag = 11 print data in file11 format
       flag = 30 print geometry to file30
disp_label is only used for geom.dat writing
-----------------------------------------------------------------------------*/

void cartesians :: print(int flag, FILE *fp_out, int new_geom_file,
char *disp_label, int disp_num) {
  int i;
  double x,y,z;
  int count = -1;

  if (flag == 0) {
    for (i = 0; i < num_atoms; ++i) {
       x = coord[++count]; y = coord[++count]; z = coord[++count];
       fprintf(fp_out,"%20.10f%20.10f%20.10f\n",x,y,z);
    }
  }
  if (flag == 1) {
     for (i = 0; i < num_atoms; ++i) {
       x = coord[++count]; y = coord[++count]; z = coord[++count];
       fprintf(fp_out,"%5.1lf%15.10f%15.10f%15.10f\n",atomic_num[i],x,y,z);
     }
  }
  if (flag == 2) {
     count = -1;
     for (i = 0; i < num_atoms; ++i) {
       x = coord[++count]; y = coord[++count]; z = coord[++count];
       fprintf(fp_out,
         "%5.1lf%15.8lf%15.10f%15.10f%15.10f\n",atomic_num[i],mass[3*i],x,y,z);
     }
  }
  if (flag == 3) {
     count = -1;
     for (i = 0; i < num_atoms; ++i) {
       x = coord[++count]; y = coord[++count]; z = coord[++count];
       fprintf(fp_out,
         "%5.1lf%15.8lf%15.10f%15.10f%15.10f\n",atomic_num[i],mass[3*i],x,y,z);
     }
     count = -1;
     for (i = 0; i < num_atoms; ++i) {
       x = grad[++count]; y = grad[++count]; z = grad[++count];
       fprintf(fp_out,
         "%35.10lf%15.10f%15.10f\n",x,y,z);
     }
  }
  if (flag == 4) {
     double *geom_out;
     if (new_geom_file) fprintf(fp_out,"%%%%\n");
     fprintf(fp_out,"%% %s\n",disp_label);
     geom_out = get_coord();
     fprintf(fp_out,"geometry%1d = (\n", disp_num);
     for (i=0;i<num_atoms;++i)
       fprintf(fp_out,"(%20.10f%20.10f%20.10f)\n",geom_out[3*i],
               geom_out[3*i+1],geom_out[3*i+2] );
     fprintf(fp_out," )\n");
     free(geom_out);
  }
  if (flag == 11) {
     fprintf(fp_out,"%s\n",disp_label);
     fprintf(fp_out,"%5d%20.10lf\n",num_atoms,energy);
     count = -1;
     for (i = 0; i < num_atoms; ++i) {
       x = coord[++count]; y = coord[++count]; z = coord[++count];
       fprintf(fp_out,
         "%20.10lf%20.10lf%20.10lf%20.10lf\n",atomic_num[i],x,y,z);
     }
     count = -1;
     for (i = 0; i < num_atoms; ++i) {
       x = grad[++count]; y = grad[++count]; z = grad[++count];
       fprintf(fp_out,"%40.10lf%20.10lf%20.10lf\n",x,y,z);
     }
  }
  if (flag == 30) {

     fprintf(outfile,"\nGeometry written to file30\n");
    fflush(outfile);
    rewind(fp_input);
    ip_set_uppercase(1);
    ip_initialize(fp_input,outfile);
    ip_cwk_clear();
    ip_cwk_add(":DEFAULT");
    ip_cwk_add(":OPTKING");
    file30_init();

     int j;
     double **geom;
     geom = init_matrix(num_atoms,3);
     for (i=0; i<num_atoms; ++i) {
        for (j=0; j<3; ++j)
           geom[i][j] = coord[3*i+j];
     }
     file30_wt_geom(geom);
     file30_close();
     ip_done();

     free_matrix(geom,num_atoms);
  }

  return;
}

