#define EXTERN
#include "includes.h"
#include "common.h"
#include <libipv1/ip_lib.h>
#include <file30.h>

char *determine_functional(void){
    int errcod;
    int depth;
    int i;
    char *exch_str,*corr_str;
    char *functional;
    
    depth = 0;
    errcod = ip_count("FUNCTIONAL",&depth,0);
    if(depth == 0){
	errcod = ip_string("FUNCTIONAL",&functional,0);
	if(errcod != IPE_OK){
	    fprintf(outfile," Must specify a functional when using ks-dft\n");
	    exit(1);
	}
    }
    else if(depth == 2){
	errcod = ip_string("FUNCTIONAL",&exch_str,1,0);
	if(errcod != IPE_OK){
	    fprintf(outfile,
		    " Exchange functional specification is invalid or missing.\n");
	    exit(1);
	}
	errcod = ip_string("FUNCTIONAL",&corr_str,1,1);
	if(errcod != IPE_OK){
	    fprintf(outfile,
		    " Correlation functional specification is invalid or missing.\n");
	    exit(1);
	}
	if (!strcmp(exch_str,"NONE"))
	    exch_str[0] = '\0';
	if (!strcmp(corr_str,"NONE"))
	    corr_str[0] = '\0';
	
	i = strlen(exch_str) + strlen(corr_str);
	if (i == 0)
	    functional = strdup("NONE");
	else {
	    functional = (char *) malloc(sizeof(char)* (i + 1));
	    sprintf(functional,"%s%s",exch_str,corr_str);
	}
	free(exch_str); free(corr_str);
    }
    else{
	fprintf(outfile,"\nwrong number of records in FUNCTIONAL keyword");
	exit(1);
    }
    return functional;
}

char *determine_grid(void){
    int errcod;
    int depth;
    int i;
    char *grid_str;
    char *r,*ang;
    char *Euler = "Euler-Mclaren / Lebedev";

    depth = 0;
    errcod = ip_count("GRID",&depth,0);
    if(depth == 0){
	grid_str = "SG-1";
	errcod = ip_string("GRID",&grid_str,0);
    }
    else if(depth == 2){
	errcod = ip_string("GRID",&r,1,0);
	errcod = ip_string("GRID",&ang,1,1);
	
	i = strlen(Euler)+strlen(r)+strlen(ang);
	
	grid_str = (char *)malloc(sizeof(char)*(i+1));
	sprintf(grid_str,"%s %s %s",Euler,r,ang);
	
	free(r);free(ang);
    }
    else{
	fprintf(outfile,"\nProblem with Grid specification: Wrong number of elements for keyword Grid");
	
	exit(1);
    }
    return grid_str;
}
