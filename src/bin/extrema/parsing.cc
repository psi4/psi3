/*###################################################################
#
#  parsing.cc
#
#  input function for extrema
#                                                     J.Kenny 7-22-00
####################################################################*/						      
#include <string.h>

extern "C" {
#include <file30.h>
#include <psio.h>
}

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







