/*#################################################################################
  coord_base.h

  coordinate base class declaration and definition
  ###############################################################################*/
#include<new>
#include<physconst.h>

double **symm_matrix_invert(double**,int,int,int);

class coord_base {
  
  protected:
    double **B, **u, **G, **A, **H, **H_old, *carts, *c_grads, *masses, *coord_arr, *grad_arr, *coord_temp,
	*coord_old, *grad_old, *coord_write;

    simple* simple_arr;

    coord_base() {
	try {simple_arr = new simple[num_simples];} 
	catch(bad_alloc) { punt("malloc error: no memory left"); }
	B = init_matrix(num_coords,3*num_atoms);
	G = init_matrix(num_coords,num_coords);
	A = init_matrix(3*num_atoms,num_coords);
	u = init_matrix(3*num_atoms,3*num_atoms);
        H = init_matrix(num_coords,num_coords);
        H_old = init_matrix(num_coords,num_coords);
        carts = init_array(3*num_atoms);
        c_grads = init_array(3*num_atoms);
        masses = init_array(num_atoms);
        coord_arr = init_array(num_coords);
        grad_arr = init_array(num_coords);
        coord_temp = init_array(num_coords);
        coord_old = init_array(num_coords);
        grad_old = init_array(3*num_atoms);
        coord_write = init_array(num_coords);

    }
    ~coord_base() { delete[] simple_arr; }
    //virtual void compute_G();

  public:
    void set_simple(int num, int type, double val, int at, int bd, int an, int tr) {
	simple_arr[num].set_simple(type, val, at, bd, an, tr); return; }
    int get_bond_atom(int i) { return simple_arr[i].get_bond(); }
    int get_angle_atom(int i) { return simple_arr[i].get_angle(); }
    int get_tors_atom(int i) {return simple_arr[i].get_tors(); }
    void set_cart(int i, double val) { carts[i]=val; return; }
    double get_cart(int i) { return carts[i]; }
    void set_c_grad(int i, double val) { c_grads[i]=val; return; }
    void set_mass(int i, double val) { masses[i]=val; return; }
    
    void print_carts() {
	int i, j;
	double **temp;
	temp = init_matrix(num_atoms,3);
	for(i=0;i<num_atoms;++i) 
	    for(j=0;j<3;++j) 
		temp[i][j] = carts[3*i+j];
	fprintf(outfile,"\nCartesian Coordinates (bohr):\n");
	print_mat(temp,num_atoms,3,outfile);
    }

    void print_c_grads() {
	int i, j;
	double **temp;
	temp = init_matrix(num_atoms,3);
	for(i=0;i<num_atoms;++i) 
	    for(j=0;j<3;++j) 
		temp[i][j] = c_grads[3*i+j];
	fprintf(outfile,"\nCartesian Gradients (hartree/bohr):\n");
	print_mat(temp,num_atoms,3,outfile);
    }
 
    virtual void compute_B() = 0;

    void print_B() {
        fprintf(outfile,"\nB matrix:\n");
        print_mat(B,num_coords,3*num_atoms,outfile);
        return;
    }

    void compute_G() {
	double **temp1;
	temp1 = init_matrix(num_coords,3*num_atoms);
	mmult(B,0,u,0,temp1,0,num_coords,3*num_atoms,3*num_atoms,0);
	mmult(temp1,0,B,1,G,0,num_coords,3*num_atoms,num_coords,0);
        free_matrix(temp1,num_coords);
	return;
    }

    void print_G() {
        fprintf(outfile,"\nG matrix:\n");
        print_mat(G,num_coords,num_coords,outfile);
        return;
    }

    void compute_A() {
	double **temp1, **temp2;
	temp1 = init_matrix(num_coords,num_coords);
	temp2 = init_matrix(3*num_atoms,num_coords);
	/* A = u Bt G-1 */
	temp1 = symm_matrix_invert(G, num_coords, 0, 1);
	mmult(B,1,temp1,0,temp2,0,3*num_atoms,num_coords,num_coords,0);
	mmult(u,0,temp2,0,A,0,3*num_atoms,3*num_atoms,num_coords,0);
        free_matrix(temp1,num_coords);
        free_matrix(temp2,3*num_atoms);
	return;
    }

    void print_A() {
        fprintf(outfile,"\nA matrix:\n");
        print_mat(A,3*num_atoms,num_coords,outfile);
        return;
    }

    void opt_step() {

	int i, j;
	double *s;
	
	s = init_array(num_coords);
  
	for(i=0;i<num_coords;++i) 
	    for(j=0;j<num_coords;++j) 
		s[i] += -H[i][j] * grad_arr[j];

	for(i=0;i<num_coords;++i) {
	    if( (fabs(s[i]) > 0.1) && (s[i] > 0.0) )
		s[i]=0.1;
	    if( (fabs(s[i]) > 0.1) && (s[i] < 0.0) )
		s[i]=-0.1;
	}
	
	fprintf(outfile,"\nNew coordinate vector:\n");
	for(i=0;i<num_coords;++i) {
            coord_write[i] = coord_arr[i];
	    coord_arr[i] += s[i];
	    fprintf(outfile,"%lf\n",coord_arr[i]);
	}
	
	free(s);
	return;
    }
 
    void print_u() {
	fprintf(outfile,"\nu matrix:\n");
        print_mat(u,3*num_atoms,3*num_atoms,outfile);
        return;
    }

    int read_opt() {

	FILE *opt_ptr;

	int i, j, error, num_elems, num_iter;

        for(i=0;i<num_coords;++i)
	    coord_temp[i] = coord_arr[i];

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
		error += ip_data("COORD","%lf",&coord_old[i],1,i);
	    }
      
	    if(error != 0)
		punt("Problem reading old coordinate values from opt.dat");
  
	    /*read old gradient vector*/
	    error += !ip_exist("GRAD",0);
	    ip_count("GRAD",&num_elems,0);
	    error += ( num_elems != num_coords );

	    for (i=0;i<num_coords;++i) {
		error += ip_data("GRAD","%lf",&grad_old[i],1,i);
	    }

	    if(error != 0)
		punt("Problem reading old gradient from opt.dat");
	    
	    /*read hmat, the inverse of the hessian*/
	    error += (!ip_exist("HMAT",0));
	    ip_count("HMAT",&num_elems,0);
	    error += ( num_elems != num_coords );
	    ip_count("HMAT",&num_elems,1,0);
	    error += ( num_elems != num_coords );
	    
	    for (i=0;i<num_coords;++i) {
		for (j=0;j<num_coords;++j) {
		    error += ip_data("HMAT","%lf", &H_old[i][j], 2, i, j);
		}
	    }
	    
	    fclose(opt_ptr);
      
	    if(error != 0)
		punt("Problem reading old hessian from opt.dat");
	}
	else num_iter = 1;

	ip_done();
	return num_iter;
}

    void grad_trans() {
	int i,j;
	for(i=0;i<num_coords;++i) {
            grad_arr[i] = 0;
	    for(j=0;j<3*num_atoms;++j) {
		grad_arr[i] += B[i][j] * c_grads[j];
	    }}}
        
    void write_opt() {

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
	    fprintf(opt_ptr,"%lf  ",grad_arr[i]);
	    ++place;
	}
	fprintf(opt_ptr,")\n\n");
	
	/*write Hessian one row at a time*/
	place=0;
	fprintf(opt_ptr,"  hmat = ( ");
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
};












