/*#################################################################################
  z_class.h

  derived z matrix class
  ###############################################################################*/

class z_class : public coord_base {

    struct z_entry* z_geom;
    double **full_geom;

  public:
    z_class();
    ~z_class() { free(z_geom); return; }
    void compute_B();
    void back_transform();
    double* cart_to_z();
    
    void update_H() {

       int i, j;
       double NUM1, NUM2, NUM12, *d_coord, *d_grad, *temp_arr1,
	   **MAT1, **temp_mat1, **temp_mat2, **MAT2, **MAT3;

       if(iteration == 1) {
	   fprintf(outfile,"\nForming empirical Hessian:\n");
	   for(i=0;i<num_coords;++i) {
	       if( simple_arr[i].get_type() == 0 )
		   H[i][i] = 1.0;
	       if( simple_arr[i].get_type() == 1 || simple_arr[i].get_type() == 2 )
		   H[i][i] = 4.0;
	   }
	   print_mat(H,num_coords,num_coords,outfile);
       }
       else {
	    
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
	       d_coord[i] = coord_arr[i] - coord_old[i];
	       d_grad[i] = grad_arr[i] - grad_old[i];
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
       
       
       
    }

    void write_file30() {
        
	double** cart_matrix;
        int i,j;
	
        cart_matrix = init_matrix(num_atoms,3);

        int pos = -1;
        for(i=0;i<num_atoms;++i) 
	    for(j=0;j<3;++j) 
		cart_matrix[i][j] = carts[++pos];

	file30_wt_geom(cart_matrix);

        pos = -1;
        for(i=1;i<num_atoms;++i) {
	    
	    if(i==1) 
		z_geom[i].bond_val = coord_arr[++pos]; 

            if(i==2) {
                z_geom[i].bond_val = coord_arr[++pos];
                z_geom[i].angle_val = coord_arr[++pos] * 180.0 / _pi;
	    }

            if(i>2) {
		z_geom[i].bond_val = coord_arr[++pos];
		z_geom[i].angle_val = coord_arr[++pos] * 180.0 / _pi;
		z_geom[i].tors_val = coord_arr[++pos] * 180.0 / _pi;
	    }
	}

        file30_wt_zmat(z_geom,num_atoms);
    }

    void opt_step() {

	int i, j;
	double *s;
	
	s = init_array(num_coords);
 
        for(i=0;i<num_coords;++i) {
	    s[i] = 0.0;
	}
  
	for(i=0;i<num_coords;++i)
	    for(j=0;j<num_coords;++j)
		s[i] += -H[i][j] * grad_arr[j];

	for(i=0;i<num_coords;++i) {
	    if( (fabs(s[i]) > 0.1) && (s[i] > 0.0) )
		s[i]=0.1;
	    if( (fabs(s[i]) > 0.1) && (s[i] < 0.0) )
		s[i]=-0.1;
	}

	for(i=0;i<num_coords;++i) {
	    fprintf(outfile,"\nequiv grp %d: %d",i, simple_arr[i].get_equiv_grp());
	    for(j=0;j<i;++j) {
		if((simple_arr[i].get_equiv_grp() == simple_arr[j].get_equiv_grp())) {
		    fprintf(outfile,"\n setting s[%d](%.20lf)=s[%d](%.20lf)",i,s[i],j,s[j]);
		    s[i] = s[j];
		}
	    }
	}

	fprintf(outfile,"\ns vector:\n");
	for(i=0;i<num_coords;++i) 
	    fprintf(outfile,"\n%.20lf",s[i]);
	
/*  	fprintf(outfile,"\nNew coordinate vector:\n"); */
/*  	for(i=0;i<num_coords;++i) { */
/*  	    coord_write[i] = coord_arr[i]; */
/*  	    coord_arr[i] += s[i]; */
/*  	    fprintf(outfile,"%lf\n",coord_arr[i]); */
/*  	} */

	fprintf(outfile,"\nNew coordinate vector:\n");
	for(i=0;i<num_coords;++i) 
            coord_write[i] = coord_arr[i];	    
	for(i=1;i<num_atoms;++i) {
	    if( (i==1) && z_geom[0].bond_opt)
		coord_arr[0] += s[0];
	    else if(i==2) {
		if(z_geom[1].bond_opt)
		    coord_arr[1] += s[1];
		//if(z_geom[1].angle_opt)
		    //coord_arr[2] += s[2];
	    }
	    else if(i>2) {
		if(z_geom[i].bond_opt) {
		    coord_arr[(i-2)*3] += s[(i-2)*3];
		    fprintf(outfile,"\ncoord %d opt",(i-2)*3);
		}
		else
		    fprintf(outfile,"\ncoord %d no opt",(i-2)*3);
		if(z_geom[i].angle_opt) {
		    coord_arr[(i-2)*3+1] += s[(i-2)*3+1];
		    fprintf(outfile,"\ncoord %d opt",(i-2)*3+1);
		}
	        else
		    fprintf(outfile,"\ncoord %d no opt",(i-2)*3+1);
		if(z_geom[i].tors_opt){
		    coord_arr[(i-2)*3+2] += s[(i-2)*3+2];
		    fprintf(outfile,"\ncoord %d opt",(i-2)*3+2);
		}
		else
		    fprintf(outfile,"\ncoord %d no opt",(i-2)*3+2);
	    }
	}

	for(i=0;i<num_coords;++i)
	    fprintf(outfile,"%.20lf\n",coord_arr[i]);
	
	free(s);
	return;
    }

};






