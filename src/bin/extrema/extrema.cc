/*###########################################################################*/
/*! \file extrema.cc
  \brief Main function and related small functions for extrema.

  Provides main function and just enough input parsing to know what
  method drivers main should call. */

  /*! \fn void main()
    \brief Main function for extrmema.
    Initializes top-level class objects and calls method driver functions. */
/*						Joseph P. Kenny 11/29/01
  ##########################################################################*/

#include "extrema.h"

int get_coord_type();



void main() {

    coord_type = get_coord_type();

    /*------------
      CARTESIANS
      ----------*/
    if(coord_type==1) {
	//carts c_obj;
    }
  

    /*---------
      ZMATRIX
      -------*/
    else if(coord_type==2) {
	zmat z_obj;
	z_obj.optimize();
    }


    /*-------------
      DELOCALIZED
      -----------*/
    else if(coord_type==3) {
	// deloc d_obj;
    }
    
    exit(0);
}



/*-----------------------------------------------------------------------------
  get_coord_type

  print intro and determine coordinate type
  ---------------------------------------------------------------------------*/

void print_intro();

int get_coord_type() {

  char *buffer;

   /*set up i/o stuff*/
  ffile(&infile,"input.dat",2);
  ffile(&outfile,"output.dat",1);

  print_intro();
 
  ip_set_uppercase(1);
  ip_initialize(infile,outfile);
  ip_cwk_clear();
  
  ip_cwk_add(":EXTREMA");
  ip_cwk_add(":DEFAULT");
  
  if(ip_exist("COORDINATES",0)) {
      errcod = ip_string("COORDINATES", &buffer,0);
      if( !strcmp(buffer,"CARTESIANS") ) {
	  fprintf(outfile,"\n  Using cartesian coordinates\n");
	  coord_type = 1; }
      else if( !strcmp(buffer,"ZMATRIX") ) {
	  fprintf(outfile,"\n  Using z-matrix coordinates\n");
	  coord_type = 2; }
      else if( !strcmp(buffer,"DELOCALIZED") ) 
	  punt("Can't do delocalized internals yet");
      else 
	  punt("Problem determining coordinate type");
      free(buffer);
  }    
  else {
      fprintf(outfile,"\n  Defaulting to z-matrix coordinates\n");
      coord_type = 2;
  }

  return coord_type;
}



void print_intro() {

     tstart(outfile);
     fprintf(outfile,"                  --------------------------------------------\n");
     fprintf(outfile,"                                   EXTREMA \n");
     fprintf(outfile,"                        Joseph P. Kenny and Rollin King \n");
     fprintf(outfile,"                  --------------------------------------------\n\n");

  return;
}








