/******************************************************************************

  zmat :: zmat

  purpose: constructor for z-matrix derived class

  parameters: none

  returns: nothing
******************************************************************************/

#include <string.h>

extern "C" {
    #include<file30.h>
}

#define EXTERN
#include "extrema.h"

zmat :: zmat()
    : internals() 
   {

  int i, j, pos;        

  try {simples = new simple[num_coords]; }
  catch(bad_alloc) { punt("error malloc'ing simples in zmat constructor"); }

  /*read z_mat from file30*/
  z_geom = file30_rd_zmat(); 
  felement = file30_rd_felement();

  /*write z_mat to the array of simple_internal*/
  for(i=1;i<num_entries;++i) {
      if( i==1 ) {
	  simples[0].set_simple(0,z_geom[1].bond_val,2,z_geom[1].bond_atom,-1,-1);
	}
      else if( i==2 ) {
	  simples[1].set_simple(0,z_geom[2].bond_val,3,z_geom[2].bond_atom,-1,-1);
	  simples[2].set_simple(1,z_geom[2].angle_val *_pi/180.0,3,z_geom[2].bond_atom,z_geom[2].angle_atom,-1);
	  j=3;
	}
      else if( i>2 ) {
	  simples[j].set_simple(0,z_geom[i].bond_val,i+1,z_geom[i].bond_atom,-1,-1);
	  simples[j+1].set_simple(1,z_geom[i].angle_val*_pi/180.0,i+1,z_geom[i].bond_atom,z_geom[i].angle_atom,-1);
	  simples[j+2].set_simple(2,z_geom[i].tors_val*_pi/180.0,i+1,z_geom[i].bond_atom,z_geom[i].angle_atom,z_geom[i].tors_atom);
	  j+=3;
	}
    }

  /*deal with symmetrically equivalent coordinates*/
  char** labels;
  labels = (char**) malloc( num_coords*sizeof(char*));
  for(i=0;i<num_coords;++i) labels[i] = (char*) malloc(20*sizeof(char));
  int p=0;
  for(i=1;i<num_entries;++i) {
      if(i==1) {
	  strcpy( labels[p], z_geom[1].bond_label );
	  ++p;
      }
      else if( i==2) {
	  strcpy(labels[p], z_geom[2].bond_label );
	  strcpy(labels[p+1], z_geom[2].angle_label );
	  p+=2;
      }
      else if (i>2) {
	  strcpy(labels[p], z_geom[i].bond_label );
	  strcpy(labels[p+1],z_geom[i].angle_label);
	  strcpy(labels[p+2],z_geom[i].tors_label);
	  p+=3;
      }
  }

  p=0;
  int is_set;
  for (i=0;i<num_coords;++i) {
      is_set=0;
      if( labels[i][0] != '\0' ) {
	  if(p==0) {
	      simples[i].set_equiv_grp(0);
	      is_set=1; ++p;
	  }
	  for(j=0;j<i;++j) {
	      if ( !strcmp(labels[i],labels[j])) {
		  is_set = 1;
		  simples[i].set_equiv_grp( simples[j].get_equiv_grp() );
	      }
	  }
	  if(!is_set) {
	      simples[i].set_equiv_grp( p );
	      ++p;
	  }
      }
      else if( labels[i][0] == '\0' )
	  simples[i].set_equiv_grp(-1);
  }

  for(i=0;i<num_coords;++i) 
      coords[i] = simples[i].get_val();

  return;
}


/***********************************************************************************
*
*  zmat :: compute_B
*
*  forms B matrix for simple internal coordinates
***********************************************************************************/
  double *B_row_bond(double* c_arr, int atom1, int atom2 );
  double *B_row_angle(double* c_arr, int atom1, int atom3, int atom2 );
  double *B_row_tors(double* c_arr, int atom1, int atom2, int atom3, int atom4 );

void zmat::compute_B() {  

    int i, j, pos=0;
    double *B_row0, *B_row1, *B_row2;
    B_row0 = init_array(3*num_entries);
    B_row1 = init_array(3*num_entries);
    B_row2 = init_array(3*num_entries);

  for(i=1;i<num_entries;++i) {
      if(i==1) {
	  B_row0 = B_row_bond(carts, i, simples[pos].get_bond()-1);
	  B[pos] = B_row0;
          ++pos;
	}
      if(i==2) {
	  B_row0 = B_row_bond(carts, i, simples[pos].get_bond()-1);
	  B_row1 = B_row_angle(carts, i, simples[pos+1].get_bond()-1, simples[pos+1].get_angle()-1);
	  B[pos] = B_row0;
	  B[pos+1] = B_row1;
	  pos += 2;
	}
      if(i>2) {
	  B_row0 = B_row_bond(carts, i, simples[pos].get_bond()-1);
	  B_row1 = B_row_angle(carts, i, simples[pos+1].get_bond()-1, simples[pos+1].get_angle()-1);
	  B_row2 = B_row_tors(carts, i, simples[pos+2].get_bond()-1, simples[pos+2].get_angle()-1, simples[pos+2].get_tors()-1);
	  B[pos] = B_row0;
	  B[pos+1] = B_row1;
	  B[pos+2] = B_row2;
	  pos += 3;
	}
    }

  /*form u*/
  for(j=0;j<num_entries;++j) {
      if(strcmp(felement[j],"X")) {
	 u[3*j][3*j] = 1.0 / masses[j]; 
	 u[3*j+1][3*j+1] = 1.0 /masses[j];
	 u[3*j+2][3*j+2] = 1.0 / masses[j];
       }
       else if (!strcmp(felement[j],"X")) {
	  fprintf(outfile,"\nMASS of X: 1.0");
	  u[3*j][3*j] = u[3*j+1][3*j+1] = u[3*j+2][3*j+2]= 1.0;
	  }
  }

  return;
}


/*###################################################################
#
#  zmat :: cart_to_internal()
#
#  computes z-mat values for cartesian coordinates
#                                                     J.Kenny 10-15-00
####################################################################*/						      
inline double norm(double* carts, int atom1, int atom2 );
inline double* unit_vec(double* carts, int atom1, int atom2 );
inline double dot_pdt( double* vec1, double* vec2 );
inline double* cross_pdt( double* vec1, double* vec2 );

void zmat :: cart_to_internal(double* z_array) {

    int i, j,pos=0;
    double *temp1, *temp2, temp_num, div, n1, n2;
    double tnum;

    temp1 = init_array(3);
    temp2 = init_array(3);



    for(i=1;i<num_entries;++i) {


	if(i==1) {
	    z_array[pos] = norm(carts, i ,simples[pos].get_bond()-1);
	    ++pos ;
	}


	if(i==2){

		z_array[pos] = norm( carts, i, simples[pos].get_bond()-1);
		temp_num = dot_pdt( unit_vec(carts, simples[pos+1].get_bond()-1, i ), 
				    unit_vec(carts, simples[pos+1].get_bond()-1, simples[pos+1].get_angle()-1 ) );
		z_array[pos+1] = acos( temp_num );
	    pos += 2;
	}


	if(i>2) {
	    
		z_array[pos] = norm( carts, i, simples[pos].get_bond()-1);
		z_array[pos+1] = acos( dot_pdt( unit_vec( carts, simples[pos+1].get_bond()-1, i ), 
					      unit_vec( carts, simples[pos+1].get_bond()-1, simples[pos+1].get_angle()-1 ) ) );

		temp1 = cross_pdt( unit_vec( carts, simples[pos+2].get_angle()-1, simples[pos+2].get_tors()-1 ), 
				   unit_vec( carts, simples[pos+2].get_angle()-1, simples[pos+2].get_bond()-1 ) );
		
		temp2 = cross_pdt( unit_vec( carts, simples[pos+2].get_bond()-1, simples[pos+2].get_angle()-1 ),
				   unit_vec( carts,simples[pos+2].get_bond()-1,i ) );

		n1 = sqrt(dot_pdt(temp1,temp1));
		n2 = sqrt(dot_pdt(temp2,temp2));

		for(j=0;j<3;++j) {
		    temp1[j] /= n1;
		    temp2[j] /= n2;
		}

		temp_num = dot_pdt(temp1,temp2);

		if(temp_num>0.999999999999999)
		    z_array[pos+2] = 0.0;
		else if(temp_num<-0.999999999999999)
		    z_array[pos+2] = _pi;
		else 
		    z_array[pos+2] = acos(temp_num);
		pos += 3;
	}
		
    }

   
    fprintf(outfile,"\nnew z-mat\n");
    for(i=0;i<num_coords;++i) {
	fprintf(outfile,"%lf\n",z_array[i]);
    }
    return;
}								     



/* inline functions */
  
inline double norm(double* carts, int atom1, int atom2 ) {

  int i;
  double norm;

  norm = 0;
  for(i=0;i<3;++i) {
      norm += ((carts[3*atom2+i] - carts[3*atom1+i]) * (carts[3*atom2+i] - carts[3*atom1+i]));
    }

  norm = sqrt(norm);

   return norm;
}


inline double* unit_vec(double* carts, int atom1, int atom2 ) {

  int i;
  double *unit_vec;

  unit_vec = init_array(3);

  for(i=0;i<3;++i) {
      unit_vec[i] = ( carts[3*atom2+i]-carts[3*atom1+i] ) / norm(carts,atom1,atom2);
    }

  return unit_vec;
}



inline double dot_pdt( double* vec1, double* vec2 ) {

  int i;
  double val;

  val=0;
  for(i=0;i<3;++i) {
      val += vec1[i] * vec2[i];
    }

  return val;
}


     
inline double* cross_pdt( double* vec1, double* vec2 ) {

  double* result_vec;

  result_vec = init_array(3);

  result_vec[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
  result_vec[1] = -(vec1[0]*vec2[2] - vec1[2]*vec2[0]);
  result_vec[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];

  return result_vec;
}


void zmat :: print_internals() {

  int i;

  /*print z-mat to output*/
  fprintf(outfile,"\n  Z-matrix:\n\n");
  for(i=-1;i<(num_entries-1);++i) {
    switch(i) {
      case -1:
	  fprintf(outfile,"%10s\n",felement[0]);
	  break;
      case 0: 
	  fprintf(outfile,"%10s %4d %14.12lf\n",
		  felement[i+1], simples[i].get_bond(), simples[i].get_val());
	  break;
      case 1: 
	  fprintf(outfile,"%10s %4d %14.12lf %4d %14.12lf\n",
		  felement[i+1], simples[i].get_bond(), simples[i].get_val(),
		  simples[i+1].get_angle(), simples[i+1].get_val()); 
	  break;
      default:
	  fprintf(outfile,"%10s %4d %14.12lf %4d %14.12lf %4d %14.12lf\n",
		  felement[i+1], simples[(i-1)*3].get_bond(), simples[(i-1)*3].get_val(),
		  simples[(i-1)*3+1].get_angle(), simples[(i-1)*3+1].get_val(),
		  simples[(i-1)*3+2].get_tors(), simples[(i-1)*3+2].get_val() );
	  break;
      }
  }
}



void zmat :: initial_H() {
    
    int i;
    
    for(i=0;i<num_coords;++i) {
	switch( simples[i].get_type() ) {
	case 0: H[i][i] = 1.0; break;
	case 1: H[i][i] = 4.0; break;
	case 2: H[i][i] = 4.0; break;
	}
    }
}
	    

void zmat :: write_file30() {

  coord_base::write_file30();
  
  int i;
  int pos = -1;
  for(i=1;i<num_entries;++i) {
      
      if(i==1) 
	  z_geom[i].bond_val = coords[++pos]; 
      
      if(i==2) {
	  z_geom[i].bond_val = coords[++pos];
	  z_geom[i].angle_val = coords[++pos] * 180.0 / _pi;
      }
      
      if(i>2) {
	  z_geom[i].bond_val = coords[++pos];
	  z_geom[i].angle_val = coords[++pos] * 180.0 / _pi;
	  z_geom[i].tors_val = coords[++pos] * 180.0 / _pi;
      }
  }
  
  file30_wt_zmat(z_geom,num_atoms);
  return;
}


void zmat :: opt_step() {
    
    int i, j;
    double *s;

    s = compute_s();

            for(i=0;i<num_coords;++i) {
            fprintf(outfile,"\nequiv grp %d: %d",i, simples[i].get_equiv_grp(
));
            for(j=0;j<i;++j) {
                if((simples[i].get_equiv_grp() == simples[j].get_equiv_grp
())) {
                    fprintf(outfile,"\n setting s[%d](%.20lf)=s[%d](%.20lf)",i,s
[i],j,s[j]);
                    s[i] = s[j];
                }
            }
        }

        fprintf(outfile,"\ns vector:\n");
        for(i=0;i<num_coords;++i)
            fprintf(outfile,"\n%.20lf",s[i]);

/*      fprintf(outfile,"\nNew coordinate vector:\n"); */
/*      for(i=0;i<num_coords;++i) { */
/*          coord_write[i] = coords[i]; */
/*          coords[i] += s[i]; */
/*          fprintf(outfile,"%lf\n",coords[i]); */
/*      } */

        fprintf(outfile,"\nNew coordinate vector:\n");
        for(i=0;i<num_coords;++i)
            coord_write[i] = coords[i];
        for(i=1;i<num_atoms;++i) {
            if( (i==1) && z_geom[0].bond_opt)
                coords[0] += s[0];
            else if(i==2) {
                if(z_geom[1].bond_opt)
                    coords[1] += s[1];
                //if(z_geom[1].angle_opt)
                    //coords[2] += s[2];
	    }
            else if(i>2) {
                if(z_geom[i].bond_opt) {
                    coords[(i-2)*3] += s[(i-2)*3];
                    fprintf(outfile,"\ncoord %d opt",(i-2)*3);
                }
                else
                    fprintf(outfile,"\ncoord %d no opt",(i-2)*3);
                if(z_geom[i].angle_opt) {
                    coords[(i-2)*3+1] += s[(i-2)*3+1];
                    fprintf(outfile,"\ncoord %d opt",(i-2)*3+1);
                }
                else
                    fprintf(outfile,"\ncoord %d no opt",(i-2)*3+1);
                if(z_geom[i].tors_opt){
                    coords[(i-2)*3+2] += s[(i-2)*3+2];
                    fprintf(outfile,"\ncoord %d opt",(i-2)*3+2);
                }
                else
                    fprintf(outfile,"\ncoord %d no opt",(i-2)*3+2);
            }
        }

        for(i=0;i<num_coords;++i)
            fprintf(outfile,"%.20lf\n",coords[i]);

        free(s);
        return;
    }



























