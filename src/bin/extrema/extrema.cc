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
void print_intro();
void start_io();
void stop_io();



int main() {

    start_io();
    print_intro();
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
	deloc d_obj;
	d_obj.optimize();
    }
    
    stop_io();
    exit(0);
    return 0;
}



/*-----------------------------------------------------------------------------
  get_coord_type

  print intro and determine coordinate type
  ---------------------------------------------------------------------------*/

int get_coord_type() {

  char *buffer;
  
  if(ip_exist("COORDINATES",0)) {
      errcod = ip_string("COORDINATES", &buffer,0);
      if( !strcmp(buffer,"CARTESIANS") ) 
	  punt("Cartesians not available");
      else if( !strcmp(buffer,"ZMATRIX") ) {
	  if( ip_exist("ZMAT",0) ) {
	      fprintf(outfile,"\n  Using z-matrix coordinates\n");
	      coord_type = ZMAT_TYPE; }
	  else
	      punt("Can't find z-matrix");
      }
      else if( !strcmp(buffer,"DELOCALIZED") ) {
	  fprintf(outfile,"\n Using delocalized internal coordinates\n");
	  coord_type = DELOC_TYPE; }
      else 
	  punt("Problem determining coordinate type");
      free(buffer);
  }    
  else {
      if( ip_exist("ZMAT",0) ) {
        fprintf(outfile,"\n  Defaulting to z-matrix coordinates\n");
        coord_type = ZMAT_TYPE; }
      else if( ip_exist("GEOMETRY",0) ) {
        fprintf(outfile,"\n  Defaulting to delocalized internal coordinates\n");
        coord_type = DELOC_TYPE; }
      else 
	punt("Problem determining coordinate type");
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



void start_io() {
    
    ffile(&infile,"input.dat",2);
    ffile(&outfile,"output.dat",1);
        
    ip_set_uppercase(1);
    ip_initialize(infile,outfile);
    ip_cwk_clear();
  
    ip_cwk_add(":EXTREMA");
    ip_cwk_add(":DEFAULT");
    ip_cwk_add(":INPUT");

    file30_init();

    return;
}



void stop_io() {

    file30_close();
    ip_done();
    tstop(outfile);
    fclose(infile);
    fclose(outfile);    
    
    return;
}
