/*###################################################################
#
#  parsing.cc
#
#  input function for extrema
#                                                     J.Kenny 7-22-00
####################################################################*/						      
#include <string.h>

#define EXTERN
#include "extrema.h"

void print_intro(void);

void parsing() {

  char *buffer;

   /*set up i/o stuff*/
  ffile(&infile,"input.dat",2);
  ffile(&outfile,"output.dat",1);

  print_intro();
 
  ip_set_uppercase(1);
  ip_initialize(infile,outfile);
  ip_cwk_clear();
  
  ip_cwk_add(":EXTREMA");
  
  if(ip_exist("COORDINATES",0)) {
      errcod = ip_string("COORDINATES", &buffer,0);
      if( !strcmp(buffer,"CARTESIANS") ) {
	  fprintf(outfile,"\n  Using cartesian coordinates\n");
	  coord_type = 1; }
      else if( !strcmp(buffer,"ZMATRIX") ) {
	  fprintf(outfile,"\n  Using z-matrix coordinates\n");
	  coord_type = 2; }
      else if( !strcmp(buffer,"DELOCALIZED") ) 
	  punt("can't do delocalized internals yet");
      else 
	  punt("CARTESIANS, ZMATRIX, or DELOCALIZED are valid values for COORDINATES");
      free(buffer);
  }    
  else {
      fprintf(outfile,"\n  Defaulting to z-matrix coordinates\n");
      coord_type = 2;
  }

  print_lvl = 1;
  errcod = ip_data("PRINT","%d",&print_lvl,0);
  grad_max = 6;
  errcod = ip_data("GRAD_MAX","%d",&grad_max,0);
  
  if(coord_type==2) {
      bond_lim = 0.1;
      errcod = ip_data("BOND_LIMIT","%lf",&bond_lim,0);
      angle_lim = 0.1;
      if(ip_exist("ANGLE_LIMIT",0)) {
	  errcod = ip_data("ANGLE_LIMIT","%lf",&angle_lim,0);
	  angle_lim = angle_lim * _pi/180.0;
      }
  }      

  return;
}


void print_intro() {

     tstart(outfile);
     fprintf(outfile,"                  --------------------------------------------\n");
     fprintf(outfile,"                             EXTREMA: an optimizer \n");
     fprintf(outfile,"                        Rollin King and Joseph P. Kenny \n");
     fprintf(outfile,"                  --------------------------------------------\n\n");

  return;
}







