extern "C" {
  #include <stdio.h>
  #include <stdlib.h>
  #include <string.h>
  #include <math.h>
  #include <libciomr.h>
  #include <ip_libv1.h>
  #include <file30.h>
  #include <physconst.h>
  #include <masses.h>
}

#define EXTERN
#include "opt.h"
#include "cartesians.h"

// This function returns the forces in cartesian coordinates in aJ/Ang
double *cartesians:: forces() {
  int i;
  double *f;
  f = init_array(3*num_atoms);

  for (i=0;i<3*num_atoms;++i)
    f[i] = -1.0 * grad[i] * _hartree2J * 1.0E18 / _bohr2angstroms;
  return f;
}

cartesians::cartesians() {
  int i, j, a, count, isotopes_given = 0, masses_given = 0;
  char label[MAX_LINELENGTH], line1[MAX_LINELENGTH];
  FILE *fp_11;
  double an,x,y,z,tval;

  if ((fp_11 = fopen("file11.dat","r")) == NULL) {
      fprintf(outfile, "Error: could not open file11.dat\n");
      exit(0) ;
  }

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

  /* read masses from input.dat or use default masses */
   count = -1;
   rewind(fp_input);
   ip_initialize(fp_input,outfile);
   ip_set_uppercase(1);
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


// flag = 0 print geometry to output.dat
// flag = 1 print geom to output.dat
// flag = 2 print geom and grad with masses
// flag = 4 print geometry to geom.dat
// flag = 30 print geometry to file30
void cartesians :: print(int flag, FILE *fp_out, int new_geom_file,
char *disp_label) {
  int i;
  double x,y,z;
  int count = -1;

//  fprintf(fp_out,"\nGeometry and Gradient read from file11.dat\n");
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
     count = -1;
     for (i = 0; i < num_atoms; ++i) {
       x = grad[++count]; y = grad[++count]; z = grad[++count];
       fprintf(fp_out,"%35.10f%15.10f%15.10f\n",x,y,z);
     }
  }
  if (flag == 4) {
     double *geom_out;
     if (new_geom_file) fprintf(fp_out,"%%%%\n");
     fprintf(fp_out,"%% %s\n",disp_label);
     geom_out = get_coord();
     fprintf(fp_out,"geometry = (\n");
     for (i=0;i<num_atoms;++i)
       fprintf(fp_out,"(%20.10f%20.10f%20.10f)\n",geom_out[3*i],
               geom_out[3*i+1],geom_out[3*i+2] );
     fprintf(fp_out," )\n");
     free(geom_out);
  }
  if (flag == 30) {

     fprintf(outfile,"\nGeometry written to file30\n");
    fflush(outfile);
    rewind(fp_input);
    ip_initialize(fp_input,outfile);
//    ip_set_uppercase(1);
//    ip_cwk_add(":DEFAULT");
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

