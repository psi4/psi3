/*#############################################################################
  coord_base.cc

  data and functions common to all coordinate systems
  ###########################################################################*/

#define EXTERN
#include"extrema.h"

extern "C" {
    #include<masses.h>
    #include<file30.h>
}

double **symm_matrix_invert(double**,int,int,int);

coord_base :: coord_base() {

    file30_init();
    num_atoms = file30_rd_natom();
    num_entries = file30_rd_nentry();
    if(coord_type==1) {
	num_entries = num_atoms;
        num_coords = 3*num_atoms;
    }
    if(coord_type==2) {
	switch(num_entries) {
	  case 2: num_coords = 1; break;
	  case 3: num_coords = 3; break;
	  default: num_coords =  (num_entries*3-6); break;
	} 
    } 
  
    coords = init_array(num_coords);
    carts = init_array(3*num_entries);
    grads = init_array(num_coords);
    c_grads = init_array(3*num_entries);
    H = init_matrix(num_coords,num_coords);
    H_old = init_matrix(num_coords,num_coords);
    u = init_matrix(3*num_entries,3*num_entries);
    coords_old = init_array(num_coords);
    grads_old = init_array(num_coords);
    coord_write = init_array(num_coords);
    coord_temp = init_array(num_coords);
    masses = init_array(num_entries);
  
    return;
}

coord_base :: ~coord_base() {
	free(coords);
	free(carts);
	free(grads);
	free(c_grads);
	free_matrix(H,num_coords);
        free_matrix(H_old,num_coords);
        free_matrix(u,3*num_entries);
	free(coords_old);
	free(grads_old);
	free(coord_write);
        free(coord_temp);
        free(masses);
	return;
    }

void coord_base :: read_carts() {

    int i; 
    double **temp;
   
    switch(coord_type) {
      case 1: 
	  temp = file30_rd_geom();
	  for(i=0;i<(3*num_atoms);++i) {
	      carts[i] = coords[i] = temp[0][i];
	  }
	  free(temp);
	  break;
      case 2: 
	temp = file30_rd_fgeom();
	for(i=0;i<(3*num_entries);++i) {
	    carts[i] = temp[0][i];
	}
	free(temp);
	break;
    }
    return;
}


void coord_base :: print_carts(double conv) {
	int i, j;
	double **temp;

	temp = init_matrix(num_entries,3);
	for(i=0;i<num_entries;++i) 
	    for(j=0;j<3;++j) 
		temp[i][j] = carts[3*i+j]*conv;
        if(conv==1.0)
	    fprintf(outfile,"\n  Cartesian Coordinates (a.u.):\n");
	else
	    fprintf(outfile,"\n  Cartesian Coordinates (angstroms):\n");
	print_mat(temp,num_entries,3,outfile);
        free_matrix(temp,num_entries);
    }

void coord_base :: print_c_grads() {
	int i, j;
	double **temp;
	temp = init_matrix(num_entries,3);
	for(i=0;i<num_entries;++i) 
	    for(j=0;j<3;++j) 
		temp[i][j] = c_grads[3*i+j];
	fprintf(outfile,"\nCartesian Gradients (hartree/bohr):\n");
	print_mat(temp,num_entries,3,outfile);
        free_matrix(temp,num_entries);
    }

void coord_base :: print_grads() {
	int i;
	switch(coord_type) {
	case 1: print_c_grads(); break;
	default: 
	    fprintf(outfile,"  \nGradient vector in internal coordinates:\n");
	for(i=0;i<num_coords;++i) 
	    fprintf(outfile,"Internal%4d: % 14.12lf\n",i,grads[i]);
	break;
	}
	return;
}
 
void coord_base :: print_u() {
	fprintf(outfile,"\nu matrix:\n");
        print_mat(u,3*num_entries,3*num_entries,outfile);
        return;
    }

void coord_base :: read_opt() {

	FILE *opt_ptr;

	int i, j, error, num_elems, num_iter;

        for(i=0;i<num_coords;++i)
	    coord_temp[i] = coords[i];

	opt_ptr = fopen("opt.dat","r");
	if( opt_ptr != NULL ) {      
	    
	    ip_set_uppercase(1);
	    ip_initialize(opt_ptr,outfile);
	    ip_cwk_add(":OPT_INFO");

	    ip_data("ITERATION","%d",&num_iter,0);
	    ++num_iter;
      
	    /*read old coordinate vector*/
	    error = 0;
	    error += !ip_exist("COORD",0);  
	    ip_count("COORD",&num_elems,0);
	    error += ( num_elems != num_coords ); 

	    for (i=0;i<num_coords;++i) {
		error += ip_data("COORD","%lf",&coords_old[i],1,i);
	    }
      
	    if(error != 0)
		punt("Problem reading old coordinate values from opt.dat");
  
	    /*read old gradient vector*/
	    error += !ip_exist("GRAD",0);
	    ip_count("GRAD",&num_elems,0);
	    error += ( num_elems != num_coords );

	    for (i=0;i<num_coords;++i) {
		error += ip_data("GRAD","%lf",&grads_old[i],1,i);
	    }

	    if(error != 0)
		punt("Problem reading old gradient from opt.dat");
	    
	    /*read hmat, the inverse of the hessian*/
	    error += (!ip_exist("HINV",0));
	    ip_count("HINV",&num_elems,0);
	    error += ( num_elems != num_coords );
	    ip_count("HINV",&num_elems,1,0);
	    error += ( num_elems != num_coords );
	    
	    for (i=0;i<num_coords;++i) {
		for (j=0;j<num_coords;++j) {
		    error += ip_data("HINV","%lf", &H_old[i][j], 2, i, j);
		}
	    }
	    
	    fclose(opt_ptr);
      
	    if(error != 0)
		punt("Problem reading old hessian from opt.dat");
	}
	else num_iter = 1;
	iteration = num_iter;

	ip_done();
}

void coord_base :: write_opt() {

	int place, i, r;
	
	FILE *opt_ptr;
	
	ffile(&opt_ptr,"opt.dat",0);
	ip_set_uppercase(1);
	ip_initialize(opt_ptr,outfile);
	ip_cwk_add("OPT_INFO");
	
	fprintf(opt_ptr,"opt_info: (\n\n");
	
	/*write iteration number*/
	fprintf(opt_ptr,"  iteration = %d\n\n",iteration);
	
	/*write coordinate vector*/
	place = 0;
	fprintf(opt_ptr,"  coord = ( ");
	for (i=0;i<num_coords;++i) {
	    if( place==8 ) {
		place = 0;
		fprintf(opt_ptr,"\n            ");
	    }
	    fprintf(opt_ptr,"%lf  ",coord_write[i]);
	    ++place;
	}
	fprintf(opt_ptr,")\n\n");
	
	/*write gradient vector*/
	place = 0;
	fprintf(opt_ptr,"  grad = ( ");
	for (i=0;i<num_coords;++i) {
	    if( place==8 ) {
		place = 0;
		fprintf(opt_ptr,"\n            ");
	    }
	    fprintf(opt_ptr,"%.20lf  ",grads[i]);
	    ++place;
	}
	fprintf(opt_ptr,")\n\n");
	
	/*write Hessian one row at a time*/
	place=0;
	fprintf(opt_ptr,"  hinv = ( ");
	for(r=0;r<num_coords;++r) {
	    if( place==0 )
		fprintf(opt_ptr,"\n         ( ");
	    for (i=0;i<num_coords;++i) {
		if( place==8 ) {
		    place = 0;
		    fprintf(opt_ptr,"\n           ");
		}
		fprintf(opt_ptr,"%lf  ",H[r][i]);
		++place;
	    }
	    fprintf(opt_ptr,")\n         ");
	    place = 0;
	}
	fprintf(opt_ptr,")\n\n");
	
	fprintf(opt_ptr,"          )\n");
	
	fclose(opt_ptr);
	ip_done();
	return;
    }


void coord_base :: opt_step() {

	int i, j;
	double *s;
	
	s = compute_s();
  
	fprintf(outfile,"\nNew coordinate vector:\n");
	for(i=0;i<num_coords;++i) {
            coord_write[i] = coords[i];
	    coords[i] += s[i];
	    fprintf(outfile,"%lf\n",coords[i]);
	}
	
	free(s);
	return;
    }

double* coord_base :: compute_s() {
    int i, j;
    double *s;
	
    s = init_array(num_coords);
    
    for(i=0;i<num_coords;++i) 
	for(j=0;j<num_coords;++j) 
	    s[i] += -H[i][j] * grads[j];
    
    for(i=0;i<num_coords;++i) {
	if( (fabs(s[i]) > 0.1) && (s[i] > 0.0) )
	    s[i]=0.1;
	if( (fabs(s[i]) > 0.1) && (s[i] < 0.0) )
	    s[i]=-0.1;
    }
    return s;
}
	

void coord_base :: update_bfgs() {

       int i, j;
       double NUM1, NUM2, NUM12, *d_coord, *d_grad, *temp_arr1,
	   **MAT1, **temp_mat1, **temp_mat2, **MAT2, **MAT3;

       /*allocate memory*/
       d_coord = init_array(num_coords);
       d_grad = init_array(num_coords);
       temp_arr1 = init_array(num_coords); 
       temp_mat1= init_matrix(num_coords,num_coords);
       temp_mat2 = init_matrix(num_coords,num_coords);
       MAT1 = init_matrix(num_coords,num_coords);
       MAT2 = init_matrix(num_coords,num_coords);
       MAT3 = init_matrix(num_coords,num_coords);
       
       /*the basic idea is to break up a nasty equation into pieces and put
	 the pieces together, it's messy*/
       
       for(i=0;i<num_coords;++i) {
	   d_coord[i] = coords[i] - coords_old[i];
	   d_grad[i] = grads[i] - grads_old[i];
       }
       
       fprintf(outfile,"\ncoordinate difference:\n");
       for(i=0;i<num_coords;++i) {
	   fprintf(outfile,"%lf\n",d_coord[i]);
       }
       
       fprintf(outfile,"\ngradient difference:\n");
       for(i=0;i<num_coords;++i) {
	   fprintf(outfile,"%lf\n",d_grad[i]);
       }
       
       for(i=0;i<num_coords;++i) {
	   for(j=0;j<num_coords;++j) {
	       temp_arr1[i] += H_old[i][j] * d_grad[j];
	   }
       }
       
       NUM1 = 0.0;
       for(i=0;i<num_coords;++i)
	   NUM1 += d_grad[i] * temp_arr1[i];
       NUM2 = 0.0;
       for(i=0;i<num_coords;++i)
	   NUM2 += d_coord[i] * d_grad[i];
       
       for(i=0;i<num_coords;++i) {
	   for(j=0;j<num_coords;++j) {
	       MAT1[i][j] = d_coord[i] * d_coord[j];
	   }
       }
       
       for(i=0;i<num_coords;++i) {
	   for(j=0;j<num_coords;++j) {
	       temp_mat1[i][j] = d_coord[i] * d_grad[j];
	   }
       }
       
       for(i=0;i<num_coords;++i) {
	   for(j=0;j<num_coords;++j) {
	       temp_mat2[i][j] = d_grad[i] * d_coord[j];
	   }
       }
       
       mmult(temp_mat1,0,H_old,0,MAT2,0,num_coords,num_coords,num_coords,0);
       
       mmult(H_old,0,temp_mat2,0,MAT3,0,num_coords,num_coords,num_coords,0);
       
       NUM12 = 1/NUM2 + NUM1/(NUM2 * NUM2);
       
       /*finally put the pieces together*/
       for(i=0;i<num_coords;++i) {
	   for(j=0;j<num_coords;++j) {
	       
	       H[i][j] = H_old[i][j] + NUM12 * MAT1[i][j]
		   - (MAT2[i][j] + MAT3[i][j]) / NUM2;
	       
	   }
       }
	   
       /*free up memory*/
       free(d_coord);
       free(d_grad);
       free(temp_arr1);
       free_matrix(temp_mat1, num_coords);
       free_matrix(temp_mat2, num_coords);
       free_matrix(MAT1, num_coords);
       free_matrix(MAT2, num_coords);
       free_matrix(MAT3, num_coords);
       
       return;
    }

void coord_base :: read_file11() {

    int i, natom, count = 1, continue_flag = 1;
    char label[MAX_LINELENGTH], line1[MAX_LINELENGTH], *tmp_ptr;
    FILE *fp_11;
    double an,x,y,z,energy;
    
    if ((fp_11 = fopen("file11.dat","r")) == NULL) {
	punt("Could not open file11.dat");
    }
    
    tmp_ptr = fgets(label, MAX_LINELENGTH, fp_11);
    
    if (tmp_ptr == NULL) {
	punt("Touble reading first line of file11.dat");
    }
    
    fgets(line1, MAX_LINELENGTH, fp_11);
    if (sscanf(line1, "%d %lf", &natom, &energy) != 2) {
	punt("Trouble reading natoms and energy from file11.dat");
    }
    
    if(natom!=num_atoms)
	punt("Number of atoms differs in file11 and file30");
    
    rewind(fp_11);
    
    while ( fgets(label, MAX_LINELENGTH, fp_11) != NULL ) {
	
	fgets(line1, MAX_LINELENGTH, fp_11);
	sscanf(line1, "%d %lf", &natom, &energy);
	
	/*read in one chunk at a time*/ 
	for (i=0; i<num_atoms ; i++) {
	    if(fscanf(fp_11, "%lf %lf %lf %lf", &an, &x, &y, &z) != 4) {
		punt("Trouble reading cartesian coordinates from file11.dat");
	    }
	    masses[i]=u[3*i][3*i]=u[3*i+1][3*i+1]=u[3*i+2][3*i+2]=an2masses[(int) an];
	    carts[3*i] = x;
	    carts[3*i+1] = y;
	    carts[3*i+2] = z;
	}
	
	for (i=0; i<num_atoms ; i++) {
	    if(fscanf(fp_11, "%lf %lf %lf", &x, &y, &z) != 3) {
		punt("Trouble reading gradients from file11.dat");
	    }
	    c_grads[3*i] = x;
	    c_grads[3*i+1] = y;
	    c_grads[3*i+2] = z;
	}
	
	fgets(line1, MAX_LINELENGTH, fp_11);
	
    }
    
    fclose(fp_11);
    return;
  }

void coord_base :: write_file30() {
    
    double** cart_matrix;
    int i,j;
    
    cart_matrix = init_matrix(num_entries,3);
    
    int pos = -1;
    for(i=0;i<num_entries;++i) 
	for(j=0;j<3;++j) 
	    cart_matrix[i][j] = carts[++pos];
    
    file30_wt_geom(cart_matrix);

    return;
}
    
















