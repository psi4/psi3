/*###########################################################################*/
/*! \file coord_base_io.cc
  \brief File read and write functions common to all coordinate types. */ 
/*						Joseph P. Kenny 11/28/01
  ###########################################################################*/

#define EXTERN
#include"extrema.h"
#include<masses.h>



/*--------------------------------------------------------------------------*/ /*! \fn coord_base::parse_input()
  \brief Parses input for info common to all coordinate types. */
/*---------------------------------------------------------------------------*/

void coord_base :: parse_input() {

    char *buffer;
    int is_set=0;

    print_lvl = NORMAL_PRINT;
    errcod = ip_data("PRINT","%d",&print_lvl,0); 
    errcod = ip_data("EXTREMA_PRINT","%d",&print_lvl,0);
    fprintf(outfile,"\n  PRINT:         %d",print_lvl);
   
    grad_max = 6;
    errcod = ip_data("GRAD_MAX","%d",&grad_max,0);
    fprintf(outfile,"\n  GRAD_MAX:      %d",grad_max);

    update = "BFGS";
    errcod = ip_string("UPDATE",&update,0);
    fprintf(outfile,"\n  UPDATE:        %s", update);
  
    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn coord_base::read_file11()
  \brief Reads gradient info from file11. */
/*---------------------------------------------------------------------------*/

void coord_base :: read_file11() {

    int i, natom, count = 1, continue_flag = 1;
    char label[133]; char line1[133]; char *tmp_ptr;
    FILE *fp_11;
    double an,x,y,z,energy, *temp_c_grads;
    temp_c_grads = init_array(3*num_atoms);
    
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
	punt("Numbers of atoms differ in file11 and file30");
    
    rewind(fp_11);
    
    while ( fgets(label, MAX_LINELENGTH, fp_11) != NULL ) {
	
	fgets(line1, MAX_LINELENGTH, fp_11);
	sscanf(line1, "%d %lf", &natom, &energy);
	
	/*read in one chunk at a time*/ 
	for (i=0; i<num_atoms ; i++) {
	    if(fscanf(fp_11, "%lf %lf %lf %lf", &an, &x, &y, &z) != 4) {
		punt("Trouble reading cartesian coordinates from file11.dat");
	    }
	    masses[i]=u[3*i][3*i]=u[3*i+1][3*i+1]=u[3*i+2][3*i+2]
		=an2masses[(int) an];

	}
	
	for (i=0; i<num_atoms ; i++) {
	    if(fscanf(fp_11, "%lf %lf %lf", &x, &y, &z) != 3) {
		punt("Trouble reading gradients from file11.dat");
	    }
	    temp_c_grads[3*i] = x;
	    temp_c_grads[3*i+1] = y;
	    temp_c_grads[3*i+2] = z;
	}
	
	fgets(line1, MAX_LINELENGTH, fp_11);
	
    }

    /* copy temp_c_grads into c_grads, putting zeroes for dummy atoms */
    int p=0;
    for(i=0;i<num_entries;++i) {
	if(!strcmp(felement[i],"X       ")) {
	    c_grads[3*i] = 0.0;
	    c_grads[3*i+1] = 0.0;
	    c_grads[3*i+2] = 0.0;
	}
	else {
	    c_grads[3*i] = temp_c_grads[3*p];
	    c_grads[3*i+1] = temp_c_grads[3*p+1];
	    c_grads[3*i+2] = temp_c_grads[3*p+2];
	    ++p;
	}
    }

    free(temp_c_grads);
    fclose(fp_11);
    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn coord_base::write_file30()
  \brief Writes new cartesians to file 30. */
/*---------------------------------------------------------------------------*/

void coord_base::write_file30() {
    
    int i,j, pos=-1;
    double** cart_matrix;
    cart_matrix = init_matrix(num_atoms,3);    

    for(i=0;i<num_entries;++i) 
	for(j=0;j<3;++j) 
	    cart_matrix[i][j] = carts[++pos];

    file30_init();
    file30_wt_geom(cart_matrix);
    file30_close();

    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn coord_base::read_opt()
  \brief Reads from opt.dat()
  
  Reads previous coordinates, gradients, and hessian inverse from opt.dat. */
/*---------------------------------------------------------------------------*/

void coord_base :: read_opt() {

	FILE *opt_ptr;

	int i, j, error, num_elems;

	opt_ptr = fopen("opt.dat","r");
	if( opt_ptr != NULL ) {      
	    
	    ip_done();
	    ip_set_uppercase(1);
	    ip_initialize(opt_ptr,outfile);
	    ip_cwk_add(":OPT_INFO");

	    ip_data("ITERATION","%d",&iteration,0);
	    ++iteration;
      
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
		    error += ip_data("HINV","%lf", &Hi_old[i][j], 2, i, j);
		}
	    }
	    
	    fclose(opt_ptr);
      
	    if(error != 0)
		punt("Problem reading old hessian from opt.dat");
	}
	else iteration=1;
	
	fprintf(outfile,"\n\n  Beginning iteration: %d\n",iteration);

	ip_done();
}



/*---------------------------------------------------------------------------*/
/*! \fn coord_base::write_opt()
  \brief Writes to opt.dat()

  Writes coordinates, gradients, and hessian inverse to opt.dat. */
/*---------------------------------------------------------------------------*/

void coord_base :: write_opt() {

    int place, i, r;
    
    FILE *opt_ptr;
    
    ffile(&opt_ptr,"opt.dat",0);
    ip_set_uppercase(1);
    ip_initialize(opt_ptr,outfile);
    
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
	    fprintf(opt_ptr,"%lf  ",Hi[r][i]);
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
