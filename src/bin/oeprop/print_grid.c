#define EXTERN
#include "includes.h"
#include "prototypes.h"
#include "globals.h"

static void print_grid_plotmtv();
static void print_grid_megapovray();
static void create_megapovray_file();
static FILE* grid_file;

void print_grid()
{
  if (grid3d)
    print_grid_megapovray();
  else
    print_grid_plotmtv();

  return;
}

void print_grid_plotmtv()
{
  int i,j,k;
  double step_x, step_y, x, y;

  fprintf(outfile," --------------------------------------------------------------\n");
  fprintf(outfile,"    *** Evaluating properties over a rectangular 2D grid ***\n");
  fprintf(outfile," --------------------------------------------------------------\n\n");
  fprintf(outfile," -Coordinates of the lower left, lower right, and upper left corners of\n");
  fprintf(outfile,"  the grid rectangle (a.u.):\n");
  fprintf(outfile,"    **            x");
  fprintf(outfile,"                     y                     z\n");
  fprintf(outfile,"   ----  --------------------");
  fprintf(outfile,"  --------------------  --------------------\n");
  fprintf(outfile,"    LL   %20.10lf  %20.10lf  %20.10lf\n",
	  grid_origin[0],grid_origin[1],grid_origin[2]);
  fprintf(outfile,"    LR   %20.10lf  %20.10lf  %20.10lf\n",
	  grid_origin[0]+grid_unit_x[0]*(grid_xy1[0]-grid_xy0[0]),
	  grid_origin[1]+grid_unit_x[1]*(grid_xy1[0]-grid_xy0[0]),
	  grid_origin[2]+grid_unit_x[2]*(grid_xy1[0]-grid_xy0[0]));
  fprintf(outfile,"    UL   %20.10lf  %20.10lf  %20.10lf\n\n\n",
	  grid_origin[0]+grid_unit_y[0]*(grid_xy1[1]-grid_xy0[1]),
	  grid_origin[1]+grid_unit_y[1]*(grid_xy1[1]-grid_xy0[1]),
	  grid_origin[2]+grid_unit_y[2]*(grid_xy1[1]-grid_xy0[1]));
  
  switch (grid) {
    case 1:
      grid_file = fopen("esp.dat","w");
      break;

    case 2:
      if (!spin_prop)
        grid_file = fopen("edens.dat","w");
      else
	grid_file = fopen("sdens.dat","w");
      break;
    
    case 3:
      if (!spin_prop)
        grid_file = fopen("edgrad.dat","w");
      else
	grid_file = fopen("sdgrad.dat","w");
      break;

    case 4:
      if (!spin_prop)
        grid_file = fopen("edlapl.dat","w");
      else
	grid_file = fopen("sdlapl.dat","w");
      break;

  }

/*    fprintf(grid_file,"%8.4lf %d %8.4lf  %8.4lf %d %8.4lf  %8.4lf %d %8.4lf\n",
            grid_xyz0[0],nix,grid_xyz1[0],
            grid_xyz0[1],niy,grid_xyz1[1],
            grid_xyz0[2],niz,grid_xyz1[2]); */
  switch (grid) {
    case 1:
    case 2:
    case 4:
      fprintf(grid_file,"$DATA = CONTOUR\n");
      fprintf(grid_file,"%% xmin = %lf xmax = %lf nx = %d\n",grid_xy0[0],grid_xy1[0],nix+1);
      fprintf(grid_file,"%% ymin = %lf ymax = %lf ny = %d\n",grid_xy0[1],grid_xy1[1],niy+1);
      fprintf(grid_file,"%% zmin = %lf zmax = %lf\n",grid_zmin,grid_zmax);
      fprintf(grid_file,"%% contfill = T meshplot = T\n");
      for(i=0;i<=niy;i++) {
        for(j=0;j<=nix;j++)
	  if (grid_pts[j][i] < grid_zmin)
            fprintf(grid_file," %lf ",grid_zmin);
	  else
	    if (grid_pts[j][i] > grid_zmax)
	      fprintf(grid_file," %lf ",grid_zmax);
	    else
	      fprintf(grid_file," %lf ",grid_pts[j][i]);
        fprintf(grid_file,"\n");
      }
      break; 
      
    case 3:
      fprintf(grid_file,"$DATA = VECTOR\n");
      fprintf(grid_file,"%% xmin = %lf xmax = %lf\n",grid_xy0[0],grid_xy1[0]);
      fprintf(grid_file,"%% ymin = %lf ymax = %lf\n",grid_xy0[1],grid_xy1[1]);
      fprintf(grid_file,"%% zmin = %lf zmax = %lf\n",0.0,0.0);
      fprintf(grid_file,"%% linecolor = 3\n");
      fprintf(grid_file,"%% xlog = off vscale = 0.2\n\n");
      step_x = (grid_xy1[0]-grid_xy0[0])/nix;
      step_y = (grid_xy1[1]-grid_xy0[1])/niy;
      for(i=0;i<=nix;i++) {
	x = grid_xy0[0] + step_x*i;
	for(j=0;j<=niy;j++) {
        y = grid_xy0[1] + step_y*j;
        if (fabs(grid_pts[i][j]) <= MAXDENSGRAD)
          fprintf(grid_file,"%9.5lf  %9.5lf  %9.5lf  %12.8lf  %12.8lf  %12.8lf\n",x,y,0.0,
                  grid_vecX[i][j],grid_vecY[i][j],0.0);
        }
      }
      break;
  }
      
  fprintf(grid_file,"$END\n");
  fclose(grid_file);
}


void print_grid_megapovray()
{
  int i,j,k;
  double step_x, step_y, step_z, x, y, z;

  fprintf(outfile," --------------------------------------------------------------\n");
  fprintf(outfile,"    *** Evaluating properties over a rectangular 3D grid ***\n");
  fprintf(outfile," --------------------------------------------------------------\n\n");
  fprintf(outfile," -Coordinates of the lower left, lower right, and upper left corners of\n");
  fprintf(outfile,"  the grid box (a.u.):\n");
  fprintf(outfile,"    **            x");
  fprintf(outfile,"                     y                     z\n");
  fprintf(outfile,"   ----  --------------------");
  fprintf(outfile,"  --------------------  --------------------\n");
  fprintf(outfile,"    LL   %20.10lf  %20.10lf  %20.10lf\n",
	  grid_origin[0],grid_origin[1],grid_origin[2]);
  fprintf(outfile,"    UR   %20.10lf  %20.10lf  %20.10lf\n",
	  grid_origin[0]+
	  grid_unit_x[0]*(grid_xyz1[0]-grid_xyz0[0])+
	  grid_unit_y[0]*(grid_xyz1[1]-grid_xyz0[1])+
	  grid_unit_z[0]*(grid_xyz1[2]-grid_xyz0[2]),
	  grid_origin[1]+
	  grid_unit_x[1]*(grid_xyz1[0]-grid_xyz0[0])+
	  grid_unit_y[1]*(grid_xyz1[1]-grid_xyz0[1])+
	  grid_unit_z[1]*(grid_xyz1[2]-grid_xyz0[2]),
	  grid_origin[2]+
	  grid_unit_x[2]*(grid_xyz1[0]-grid_xyz0[0])+
          grid_unit_y[2]*(grid_xyz1[1]-grid_xyz0[1])+
          grid_unit_z[2]*(grid_xyz1[2]-grid_xyz0[2]));

  /*--- Write out a data file ---*/
  switch (grid) {
  case 5:
  case 6:
      grid_file = fopen("mo.dat","w");
      break;
  }

  switch (grid) {
  case 5:
  case 6:
    for(k=0;k<=niz;k++)
      for(j=0;j<=niy;j++)
	for(i=0;i<=nix;i++) {
	  fprintf(grid_file,"%15.8lf\n",grid3d_pts[i][j][k]);
	}
      break;
  }

  fclose(grid_file);

  /*--- Write out a command file ---*/
  create_megapovray_file();
  
  return;
}


void create_megapovray_file()
{
  int i, j, z;
  double radius, midx, midy, midz;
  double camera_distance, frame_width, frame_height, frame_center[3];
  double dimx, dimy, dimz, maxdim;
  FILE *mpvfile;
#include <rgb.h>

  dimx = grid_xyz1[0] - grid_xyz0[0];
  dimy = grid_xyz1[1] - grid_xyz0[1];
  dimz = grid_xyz1[2] - grid_xyz0[2];
#define MAX(x,y) (((x) > (y)) ? (x) : (y))
  maxdim = MAX(dimz,MAX(dimx,dimy));
  frame_width = maxdim + 2.0;
  frame_height = maxdim + 2.0;
  frame_center[0] = grid_origin[0] + 0.5*dimx;
  frame_center[1] = grid_origin[1] + 0.5*dimy;
  frame_center[2] = grid_origin[2] + 0.5*dimz;
  camera_distance = 3.0*frame_height;
    
  compute_connectivity();

  mpvfile = fopen("mo.pov","w");
  fprintf(mpvfile,"// File: mo.pov\n");
  fprintf(mpvfile,"// Creator: oeprop (Psi 3.0)\n");
  fprintf(mpvfile,"// Version: MegaPov 0.5\n\n");

  fprintf(mpvfile,"#version unofficial MegaPov 0.5\n\n");

  fprintf(mpvfile,"// Camera\n");
  fprintf(mpvfile,"camera {\n");
  fprintf(mpvfile,"    orthographic\n");
  fprintf(mpvfile,"    location <%lf, %lf, %lf>\n",
	  camera_distance,frame_center[1],frame_center[2]);  
  fprintf(mpvfile,"    look_at <%lf, %lf, %lf>\n",
	  frame_center[0],frame_center[1],frame_center[2]);
  fprintf(mpvfile,"    up <%lf, %lf, %lf>\n",
	  frame_center[0],frame_center[1],frame_center[2]+frame_height);
  fprintf(mpvfile,"    right <%lf, %lf, %lf>\n",
	  frame_center[0],frame_center[1]+frame_width,frame_center[2]);
  fprintf(mpvfile,"}\n\n");

  fprintf(mpvfile,"// Light\n");
  fprintf(mpvfile,"light_source {<%lf, %lf, %lf> color rgb <1.0, 1.0, 1.0>}\n\n",
	  camera_distance*10.0,camera_distance*(-10.0),camera_distance*(-10.0));

  fprintf(mpvfile,"// Background\n");
  fprintf(mpvfile,"background { color rgb <1.0, 1.0, 1.0>}\n\n");

  fprintf(mpvfile,"union {\n");
  fprintf(mpvfile,"// Objects\n");
  fprintf(mpvfile,"  // union finish\n");
  fprintf(mpvfile,"  #declare F = finish {specular 0.4 roughness 0.005 diffuse 0.8 ambient 0.2}\n");
  fprintf(mpvfile,"  // transparency of atomic surfaces \n");
  fprintf(mpvfile,"  #declare T = 0\n");
  fprintf(mpvfile,"  // no_shadow\n");
  fprintf(mpvfile,"  // hollow\n\n");
  
  fprintf(mpvfile,"  // Atoms\n");
  for(i=0;i<natom;i++) {
    z = (int)zvals[i];
    fprintf(mpvfile,"  // Atom # %d, charge %d\n",i,z);
    fprintf(mpvfile,"  object {\n");
    fprintf(mpvfile,"    sphere { < %lf, %lf, %lf > ",geom[i][0],geom[i][1],geom[i][2]);
    if (z <= 2)
      radius = 0.4;
    else if (z <= 10)
      radius = 0.5;
    else
      radius = 0.6;
    fprintf(mpvfile,"%lf }\n",radius);
    fprintf(mpvfile,"     texture {\n");
    if (z <= LAST_RGB_INDEX)
      fprintf(mpvfile,"       pigment {color rgbt <%lf, %lf, %lf, T>}\n",atomic_rgb[z][0],atomic_rgb[z][1],atomic_rgb[z][2]);
    else
      fprintf(mpvfile,"       pigment {color rgbt <%lf, %lf, %lf, T>}\n",atomic_rgb[0][0],atomic_rgb[0][1],atomic_rgb[0][2]);
    fprintf(mpvfile,"       finish {F}\n");
    fprintf(mpvfile,"     }\n");
    fprintf(mpvfile,"  }\n");
  }
  fprintf(mpvfile,"\n");

  fprintf(mpvfile,"  // Bonds\n");
  for(i=0;i<natom;i++)
    for(j=0;j<i;j++)
      if (connectivity[i][j]) {
	fprintf(mpvfile,"  // Bond between atoms %d and %d\n",i,j);
	midx = 0.5*(geom[i][0] + geom[j][0]);
	midy = 0.5*(geom[i][1] + geom[j][1]);
	midz = 0.5*(geom[i][2] + geom[j][2]);
	fprintf(mpvfile,"  object {\n");
	fprintf(mpvfile,"    cylinder { < %lf, %lf, %lf > < %lf, %lf, %lf > 0.15 }\n",
		geom[i][0],geom[i][1],geom[i][2],
		midx,midy,midz);
	fprintf(mpvfile,"     texture {\n");
	z = (int) zvals[i];
	if (z <= LAST_RGB_INDEX)
	    fprintf(mpvfile,"       pigment {color rgbt <%lf, %lf, %lf, T>}\n",atomic_rgb[z][0],atomic_rgb[z][1],atomic_rgb[z][2]);
	else
	    fprintf(mpvfile,"       pigment {color rgbt <%lf, %lf, %lf, T>}\n",atomic_rgb[0][0],atomic_rgb[0][1],atomic_rgb[0][2]);
	fprintf(mpvfile,"       finish {F}\n");
	fprintf(mpvfile,"     }\n");
	fprintf(mpvfile,"  }\n");
	fprintf(mpvfile,"  object {\n");
	fprintf(mpvfile,"    cylinder { < %lf, %lf, %lf > < %lf, %lf, %lf > 0.15 }\n",
		geom[j][0],geom[j][1],geom[j][2],
		midx,midy,midz);
	fprintf(mpvfile,"     texture {\n");
	z = (int) zvals[j];
	if (z <= LAST_RGB_INDEX)
	    fprintf(mpvfile,"       pigment {color rgbt <%lf, %lf, %lf, T>}\n",atomic_rgb[z][0],atomic_rgb[z][1],atomic_rgb[z][2]);
	else
	    fprintf(mpvfile,"       pigment {color rgbt <%lf, %lf, %lf, T>}\n",atomic_rgb[0][0],atomic_rgb[0][1],atomic_rgb[0][2]);
	fprintf(mpvfile,"       finish {F}\n");
	fprintf(mpvfile,"     }\n");
	fprintf(mpvfile,"  }\n");
      }
  
  fprintf(mpvfile,"  #declare ISO = 1;\n");
  fprintf(mpvfile,"  #declare NX = %d;\n",nix);
  fprintf(mpvfile,"  #declare NY = %d;\n",niy);
  fprintf(mpvfile,"  #declare NZ = %d;\n",niz);
  fprintf(mpvfile,"  #declare RX = %lf;\n",grid_origin[0]);
  fprintf(mpvfile,"  #declare RY = %lf;\n",grid_origin[1]);
  fprintf(mpvfile,"  #declare RZ = %lf;\n",grid_origin[2]);
  fprintf(mpvfile,"  #declare LX = %lf;\n",grid_xyz1[0] - grid_xyz0[0]);
  fprintf(mpvfile,"  #declare LY = %lf;\n",grid_xyz1[1] - grid_xyz0[1]);
  fprintf(mpvfile,"  #declare LZ = %lf;\n",grid_xyz1[2] - grid_xyz0[2]);
  fprintf(mpvfile,"  #declare LEVSCALE = 11;\n\n");

  fprintf(mpvfile,"  #ifdef (ISO)\n");
  fprintf(mpvfile,"  #declare wfun1 = function {\"data_3D_3\",\n");
  fprintf(mpvfile,"     <-LEVSCALE>, library \"i_dat3d\",  \"mo.dat\", <NX+1,NY+1,NZ+1,0>}\n");
  fprintf(mpvfile,"  #declare wfun2 = function {\"data_3D_3\",\n");
  fprintf(mpvfile,"      <LEVSCALE>, library \"i_dat3d\",  \"mo.dat\", <NX+1,NY+1,NZ+1,0>}\n");
  fprintf(mpvfile,"  union {\n");
  fprintf(mpvfile,"    isosurface {\n");
  fprintf(mpvfile,"      function {wfun1}\n");
  fprintf(mpvfile,"      contained_by { box {<0,0,0>,<NX,NY,NZ>} }\n");
  fprintf(mpvfile,"      threshold -1\n");
  fprintf(mpvfile,"      max_gradient 1.000\n");
  fprintf(mpvfile,"      eval\n");
  fprintf(mpvfile,"      texture { pigment { color rgbt <1.0,0.0,0.0,0.8> } }\n");
  fprintf(mpvfile,"    }\n");
  fprintf(mpvfile,"    isosurface {\n");
  fprintf(mpvfile,"      function {wfun2}\n");
  fprintf(mpvfile,"      contained_by { box {<0,0,0>,<NX,NY,NZ>} }\n");
  fprintf(mpvfile,"      threshold -1\n");
  fprintf(mpvfile,"      max_gradient 1.000\n");
  fprintf(mpvfile,"      eval\n");
  fprintf(mpvfile,"      texture { pigment { color rgbt <0.0,1.0,0.0,0.8> } }\n");
  fprintf(mpvfile,"    }\n");
  fprintf(mpvfile,"    translate <RX*NX/LX,RY*NY/LY,RZ*NZ/LZ>\n");
  fprintf(mpvfile,"    scale <LX/NX,LY/NY,LZ/NZ>\n");
  fprintf(mpvfile,"  }\n");
  fprintf(mpvfile,"  #end\n\n");
  
  fprintf(mpvfile,"  // Rotate\n");
  fprintf(mpvfile,"  rotate 0*z\n");
  fprintf(mpvfile,"}\n");

  fclose(mpvfile);

  return;
}
