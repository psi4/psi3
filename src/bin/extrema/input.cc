/*###################################################################
#
#  input.cc
#
#  input function for extrema
#                                                     J.Kenny 7-22-00
####################################################################*/						      

#include <stdio.h>
#include <stdlib.h>

extern "C" {
#include <libciomr.h>
#include <ip_libv1.h>
#include <file30.h>
#include <physconst.h>
#include <psio.h>
}

#define EXTERN
#include "opt.h"

void input() {

  int i, j;  
  
  /*set up i/o stuff*/
  ffile(&infile,"input.dat",2);
  ffile(&outfile,"output.dat",1);
 
  ip_set_uppercase(1);
  ip_initialize(infile,outfile);
  ip_cwk_clear();
  
  /*read stuff from file30*/
  file30_init();
  num_atoms = file30_rd_natom();
  file30_close;

  ip_cwk_add(":INPUT");
  if(ip_exist("ZMAT",0)) { 
     coord_type=0;
  }
  else punt("z-matrix not found in input ... can't do whatever it is you want to do yet");
   
  if(coord_type==0) {
      if( num_atoms==2 )
	  num_simples = 1;
      else if( num_atoms>2 )
	  num_simples = 3*num_atoms - 6;
      else punt("Bad number of atoms or trouble reading from input.dat");
      num_coords = num_simples;
  }
}
  





