/* This file contains functions such as Rotation, distance, and bond angle
calculations */

#define EXTERN
#include <stdio.h>
#include <libciomr.h>
#include <stdlib.h>
#include <math.h>
#include "input.h"
#include <physconst.h>
#include "global.h"
#include "defines.h"
#include <qt.h>

/* Functions */

/*
   strange, but necessary, gprgid function
*/
char *gprgid()
{
   char *prgid = "INPUT";
   return(prgid);
}

/*
   quit routine
*/
void punt(const char *msg)
{
  printf(msg);
  exit(1);
}


/* 
   UNIT_VECTORS(): Function calculates unit vectors between each
   pair of atoms and stores them in E. 
*/
void unit_vectors(int num_atoms, double *X, double *Y, double *Z, 
                  double **Distance, double ***E )
{
   int i, j ;

   for(i=0; i<num_atoms; i++) {
      for (j=0; j<num_atoms; j++) {
         if (i != j) {
            E[i][j][0] = -(X[i] - X[j]) / (Distance[i][j]);
            E[i][j][1] = -(Y[i] - Y[j]) / (Distance[i][j]);
            E[i][j][2] = -(Z[i] - Z[j]) / (Distance[i][j]);
            }
         else {
            E[i][j][0] = 0.0 ;
            E[i][j][1] = 0.0 ;
            E[i][j][2] = 0.0 ;
           }
        }
     }
}


/* 
   UNIT_VEC(): Function calculates the unit vector between a
   pair of atoms and stores it in AB. 
*/

void unit_vec(double *B, double *A, double *AB)
{  double norm = 0.0;
   int i;

   for(i=0;i<3;i++)
     norm += (A[i]-B[i])*(A[i]-B[i]);
   norm = sqrt(norm);
   for(i=0;i<3;i++)
     AB[i] = (B[i] - A[i]) / norm;
   return;
}        


/*
** CALC_BOND_ANGLES(): Function calculates all non-redundant bond angles
**   between all combinations of three atoms, and writes to file fpo.
**   Note: the first index of BondAngles is always smaller; i.e.
**   angle 5-3-1 is stored only as 1-3-5
*/
void calc_bond_angles(double ***E, 
                      double ***Bond_Angle, 
                      double **Distance, FILE *outfile)
{
int i, j, k ;
double dotprod ;
double angle ;

   fprintf(outfile, "\nBond Angles:\n") ;
   for (i=0; i<num_atoms; i++) {
      for (j=0; j<num_atoms; j++) {
         for (k=i; k<num_atoms; k++) {
            if ( (i != j) && (i != k) && (j != k) ) {
               dotprod = dot_prod(E[j][i], E[j][k]) ;
               if (dotprod > 1.00000) angle = 0.0000 ;
               else if (dotprod < -1.00000) angle = _pi ;
               else angle = acos(dotprod) ;
               Bond_Angle[i][j][k] = angle ;
               if (((Distance[i][j] < 4.000) && (Distance[j][k] < 4.000))){ 
                  fprintf(outfile, "%2d-%2d-%2d    %14.8lf\n", i+1, j+1, k+1,
                     (angle*180.00/_pi));
                  }
               }
            }
         }
      }
}

/*
   Function to calculate all internuclear distances 
*/
void calc_distance(double **geom, double *A, int num)
{

   int i = 0;
   int j = 0;
   double dist;
   double temp1=0.0;
   double temp2=0.0;
   double temp3=0.0;

   for(i=1;i<num;i++){
      for(j=0;j<i;j++){
         temp1  = (geom[i][0]-geom[j][0])*(geom[i][0]-geom[j][0]); 
         temp2  = (geom[i][1]-geom[j][1])*(geom[i][1]-geom[j][1]); 
         temp3  = (geom[i][2]-geom[j][2])*(geom[i][2]-geom[j][2]); 
         dist = sqrt(temp1 + temp2 + temp3);
	 if (dist < ZERO_BOND_DISTANCE) {
	   printf("\nERROR: Atoms %d and %d are separated by only %lf!\n\n",i+1,j+1,dist);
	   exit(1);
	 }
	 else
	   A[ioff[i]+j] = dist;
        }
     } 
}


/*
   Function to calculate the nuclear repulsion energy
*/
void Nuc_repulsion(double *Distance, double *repulsion)
{

   int i, j;
   *repulsion = 0.0;

   for(i=1;i<num_atoms;i++)
      for(j=0;j<i;j++)
         *repulsion += nuclear_charges[i]*nuclear_charges[j]/Distance[ioff[i]+j] ;
   return;
}

/*
   Function to take the dot product
*/

double dot_prod(double *v1, double *v2)
{
   return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}


/*
   Function to compute the cross product of 2 vectors
*/

void cross_prod(double *v1, double *v2, double *out)
{
   out[0] =     v1[1]*v2[2]-v1[2]*v2[1];
   out[1] =    -v1[0]*v2[2]+v1[2]*v2[0];
   out[2] =     v1[0]*v2[1]-v1[1]*v2[0];
   return;
}


/*
  Function returns 1 if vectors are equivalent and 0 if they are not
*/

int vect_equiv(double *A, double *B)
{
  if ((fabs(A[0] - B[0]) < ZERO) && (fabs(A[1] - B[1]) < ZERO) && (fabs(A[2] - B[2]) < ZERO))
    return 1;
  else
    return 0;
}


/*
  Function returns unit vector - median of vectors A and B
*/

void median_vec(double *A, double *B, double *median)
{
  double norm, x, y, z;

  x = (A[0] + B[0])/2;
  y = (A[1] + B[1])/2;
  z = (A[2] + B[2])/2;
  
  norm = sqrt(x*x + y*y + z*z);
  median[0] = x/norm;
  median[1] = y/norm;
  median[2] = z/norm;
  return;
}


/*
  Function returns 1 if vectors are inversion symmetry related, 0 otherwise
*/

int inv_related(double *A, double *B)
{
  if ((fabs(A[0] + B[0]) < ZERO) && (fabs(A[1] + B[1]) < ZERO) && (fabs(A[2] + B[2]) < ZERO))
    return 1;
  else
    return 0;
}


/*
  Function rotates geometry[][] by new_coord[][]
*/

void rotate_geometry(double **geom, double **new_coord)
{
  double **new_geom;
  int i;

  new_geom = init_matrix(num_atoms, 3);

  mmult(geom,0,new_coord,0,new_geom,0,num_atoms,3,3,0);
  for(i=0;i<num_atoms;i++) {
    geom[i][0] = new_geom[i][0];
    geom[i][1] = new_geom[i][1];
    geom[i][2] = new_geom[i][2];
  }

  free_matrix(new_geom, num_atoms);
  return;
}


/*
  Function forms a (3x3) matrix out of 3 3-vectors (aligned in columns)
*/

void vectors_to_matrix(double *v1, double *v2, double *v3, double **matrix)
{
  matrix[0][0] = v1[0]; matrix[1][0] = v1[1]; matrix[2][0] = v1[2];
  matrix[0][1] = v2[0]; matrix[1][1] = v2[1]; matrix[2][1] = v2[2];
  matrix[0][2] = v3[0]; matrix[1][2] = v3[1]; matrix[2][2] = v3[2];
  return;
}


/*Functions to allocate memory for character arrays*/

char *init_char_array(int B)
{
   char *A;
   int i;
   A = (char *)malloc((B+1)*sizeof(char));
   for(i=0;i<B;i++){
    A[i] = ' ';
   }
   A[B] = '\0';
  return A;
}

char **init_char_matrix(int A, int B)
{
   char **mat;
   int i,j;

   mat = (char **)malloc(A*sizeof(char *));

   for(i = 0;i<A;i++){
      mat[i] = (char *)malloc(B*sizeof(char));
     }
   for(i=0;i<A;i++)
     mat[i] = init_char_array(B);
   return mat;
}

void free_char_matrix(char **mat, int A)
{
  int i;

  for(i=0;i<A;i++)
    free(mat[i]);
  free(mat);
}


void setup(maxdf)
  int maxdf;
{
  int i;

  if(maxdf<8) maxdf = 8;
  df = (int *) malloc(sizeof(int)*maxdf);
  ioff = (int *) malloc(sizeof(int)*maxdf);

/* df[i] = (i-1)!! */
  df[0] = 1;
  df[1] = 1;
  df[2] = 1;
  for(i=3; i<maxdf; i++)
    df[i] = (i-1)*df[i-2];

  ioff[0] = 0;
  for(i=1;i<maxdf;i++)
    ioff[i] = ioff[i-1]+i;
}


double int_pow(a, p)
  double a;
  int p;
{
  register int i;
  double b = 1.0;

  for(i=0; i<p; i++) b = b*a;
  return b;
}
