#define EXTERN
#include <stdio.h>
#include <libciomr/libciomr.h>
#include <stdlib.h>
#include <math.h>
#include "input.h"
#include "global.h"
#include "defines.h"
#include <physconst.h>
#include <masses.h>

void reorient()
{
  int i,j;
  int degen;
  int deg_IM1, deg_IM2;
  int nspher_set;
  int prox_i, prox_j;
  double Xcm = 0.0;
  double Ycm = 0.0;
  double Zcm = 0.0;
  double mass = 0.0;
  double tmp, abs, rel;
  double Zmax, r, min_ij;
  double **IT, *IM, **ITAxes;
  double *v1, *v2, *v3, **R;
  double **sset_geom, *sset_dist;
  double origin[] = {0.0, 0.0, 0.0};

  const double im2rotconst = 0.25/(M_PI*_c_au*_bohr2cm);

    v1 = init_array(3);
    v2 = init_array(3);
    v3 = init_array(3);
    IT = block_matrix(3,3);
    IM = init_array(3);
    ITAxes = block_matrix(3,3);

    for(i=0;i<num_atoms;i++) {
      tmp = an2masses[(int)nuclear_charges[i]];
      Xcm += tmp*geometry[i][0];
      Ycm += tmp*geometry[i][1];
      Zcm += tmp*geometry[i][2];
      mass += tmp;
    }
    Xcm /= mass; Ycm /= mass; Zcm /= mass;

    if (!no_comshift) {
      /*full geom center-of-mass shift*/
      for(i=0;i<num_allatoms;i++) {
	full_geom[i][0] -= Xcm;
	full_geom[i][1] -= Ycm;
	full_geom[i][2] -= Zcm;
      }
    }


    /*Computing inertia tensor, moments of inertia, and principal axes*/
    if (num_atoms > 1) {
      for(i=0;i<num_atoms;i++) {
        tmp = an2masses[(int)nuclear_charges[i]]/_au2amu;
        IT[0][0] += tmp*(geometry[i][1]*geometry[i][1] + geometry[i][2]*geometry[i][2]);
        IT[1][1] += tmp*(geometry[i][0]*geometry[i][0] + geometry[i][2]*geometry[i][2]);
        IT[2][2] += tmp*(geometry[i][0]*geometry[i][0] + geometry[i][1]*geometry[i][1]);
	IT[0][1] -= tmp*geometry[i][0]*geometry[i][1];
	IT[0][2] -= tmp*geometry[i][0]*geometry[i][2];
	IT[1][2] -= tmp*geometry[i][1]*geometry[i][2];
      }
      IT[1][0] = IT[0][1];
      IT[2][0] = IT[0][2];
      IT[2][1] = IT[1][2];
      sq_rsp(3,3,IT,IM,1,ITAxes,1.0E-14);
      IM[0] = fabs(IM[0]); /*Fixing numerical errors in the linear case*/
      fprintf(outfile,"\n  -Rotational constants (cm-1) :\n");
      if (IM[0] < ZERO_MOMENT_INERTIA) /* Linear molecule */
	fprintf(outfile,"    A = **********  ");
      else   /* Regular molecule */
	fprintf(outfile,"    A = %10.5lf  ",im2rotconst/IM[0]);
      if (IM[1] < ZERO_MOMENT_INERTIA)  /* Atom */
	fprintf(outfile,"B = **********  C = **********\n");
      else /* molecule */
	fprintf(outfile,"B = %10.5lf  C = %10.5lf\n",im2rotconst/IM[1],im2rotconst/IM[2]);

    
      /*Ensuring the righthandedness of the reference coordinate system*/
      v1[0] = ITAxes[0][1];
      v1[1] = ITAxes[1][1];
      v1[2] = ITAxes[2][1];
      v2[0] = ITAxes[0][2];
      v2[1] = ITAxes[1][2];
      v2[2] = ITAxes[2][2];
      cross_prod(v1,v2,v3);
      ITAxes[0][0] = v3[0];
      ITAxes[1][0] = v3[1];
      ITAxes[2][0] = v3[2];

      /*Computing degeneracy*/
      degen = 0;
      for(i=0;i<2;i++)
	for(j=i+1;j<3 && degen<2;j++) {
	  abs = fabs(IM[i] - IM[j]);
	  tmp = (IM[i] > IM[j]) ? IM[i] : IM[j];
	  if (abs > 1.0E-14)
            rel = abs/tmp;
          else
	    rel = 0.0;
	  if (rel < ZERO_MOMENT_INERTIA) {
	    degen++;
	    deg_IM1 = i;
	    deg_IM2 = j;
	  }
	}

      /*Reorient if degen < 2 (non-spherical top case).
	Otherwise hope user knows what he/she's doing and leave it as it is.*/
      if (degen < 2 && !no_reorient) {
	rotate_full_geom(ITAxes);
      }

      /*If degen=0 (asymmetric top) - do nothing
	   degen=1 (linear molecule or symmetric top) - set unique axis along Z
           denen=2 (atom or spherical top) - do lots of stuff - see below*/
      if (degen == 0) {
	fprintf(outfile,"    It is an asymmetric top.\n");
	rotor = asymmtop;
      }
      else if (degen == 1)
	switch (deg_IM1 + deg_IM2) {
	  case 3: /*B and C are degenerate - linear or prolate symm. top.
		    A is the unique axis. Rotate around y-axis.*/
                  if (!no_reorient) {
	            R = block_matrix(3,3);
		    R[1][1] = 1.0;
		    R[2][0] = -1.0;
		    R[0][2] = 1.0;
		    rotate_full_geom(R);
		    free_block(R);
                  }
		  if (IM[0] < ZERO_MOMENT_INERTIA) {
		    fprintf(outfile,"    It is a linear molecule.\n");
		    rotor = linear;
		  }
		  else {
		    fprintf(outfile,"    It is a prolate symmetric top.\n");
		    rotor = symmtop;
		  }
		  break;
	  case 1: /*A and B are degenerate - oblate top.
		    C is the unique axis. Do nothing.*/
	          fprintf(outfile,"    It is an oblate symmetric top.\n");
		  rotor = symmtop;
	          break;
	}
      else if (degen == 2) { /*Man, this piece of code is nasty!!!*/
	fprintf(outfile,"    It is a spherical top.\n");
	rotor = sphtop;
	Zmax = 0;
	for(i=0;i<num_atoms;i++) /*Finding the heaviest type of atoms positioned NOT in the origin*/
	  if (Zmax < (int)nuclear_charges[i] && sqrt(dot_prod(geometry[i],geometry[i])) > ZERO)
	    Zmax = (int)nuclear_charges[i];
	/*Find subset of heaviest atoms at the same distance from the origin - so called spherical set
	  and store their geometries in sset_geom[][] */
	r = 0.0;
	for(i=0;i<num_atoms && r < ZERO;i++)
	  r = (nuclear_charges[i] == Zmax) ? sqrt(dot_prod(geometry[i],geometry[i])) : 0.0;
	sset_geom = init_matrix(num_atoms,3);
	sset_geom[0][0] = geometry[i-1][0];
	sset_geom[0][1] = geometry[i-1][1];
	sset_geom[0][2] = geometry[i-1][2];
	nspher_set = 1;
	for(j=i;j<num_atoms;j++)
	  if (nuclear_charges[j] == Zmax) {
	    tmp = sqrt(dot_prod(geometry[i],geometry[i]));
	    if ( fabs(tmp - r) < ZERO) {
	      sset_geom[nspher_set][0] = geometry[j][0];
	      sset_geom[nspher_set][1] = geometry[j][1];
	      sset_geom[nspher_set][2] = geometry[j][2];
	      nspher_set++;
	    }
	  }
	/*Find the unique orientation. If nspher_set = 4, 6, or 8 - we are dealing with a tetrahedron,
	  octahedron or cube respectively - and treat these special cases rather easily.
	  In general case, the procedure is as follows:
	  1. Pick an atom (1)
	  2. Pick an atom (2) not related to (1) by inversion
	  3. Loop over the spherical set and find all atoms {(3)}
	     so that dist(1)-(2) == dist(1)-(3)
	  4. switch(1+Number of atoms in {(3)})
	       1,2 - pick different (2)
	       3 - atom(1) is on a C3 axis and (2) and {(3)} are related by C3
	       4 - either atom(1) is on a C2 axis or there're two C3's next to each other
	       5 - find (4) in {(3)} that dist(1)-(4) == dist(1)-(2)
	           1, 2, 4 are related and perpendicular to C3
	  5. Do that again trying to find either 2 C3 axes or C3 and C2 or two C2's
	  6. Let's go ... */
        if (!no_reorient)
	switch(nspher_set) {
	  case 4: /*Tetrahedron*/
	          median_vec(sset_geom[0], sset_geom[1], v1);
		  median_vec(sset_geom[1], sset_geom[2], v2);
		  cross_prod(v1, v2, v3);
		  vectors_to_matrix(v1, v2, v3, ITAxes);
		  rotate_full_geom(ITAxes);
		  break;

	  case 6: /*Octahedron*/
  	          if (!inv_related(sset_geom[0], sset_geom[1]))
		    median_vec(sset_geom[0], sset_geom[1], v1);
		  else 
		    median_vec(sset_geom[0], sset_geom[2], v1);
		  unit_vec(sset_geom[0],origin,v2);
		  cross_prod(v1, v2, v3);
		  unit_vec(v3,origin,v3);
		  cross_prod(v3, v1, v2);
		  vectors_to_matrix(v1, v2, v3, ITAxes);
		  rotate_full_geom(ITAxes);
	          break;

	  case 8: /*Cube*/
	          sset_dist = init_array(nspher_set*(nspher_set+1)/2);
		  calc_distance(sset_geom,sset_dist,nspher_set);
	          min_ij = sset_dist[ioff[1]+0];
		  prox_i = 1;
		  for(i=2;i<nspher_set-2;i++)
		    if (min_ij > (tmp = sset_dist[ioff[i]+0]) && !inv_related(sset_geom[0],sset_geom[i])) {
		      min_ij = tmp;
		      prox_i = i;
		      break;
		    }
		  for(j=prox_i+1;j<nspher_set;j++)
		    if (fabs(min_ij - sset_dist[ioff[j]+0]) < ZERO) {
		      prox_j = j;
		      break;
		    }
		  unit_vec(sset_geom[0], sset_geom[prox_i], v1);
		  unit_vec(sset_geom[0], sset_geom[prox_j], v2);
		  cross_prod(v1, v2, v3);
		  vectors_to_matrix(v1, v2, v3, ITAxes);
		  rotate_full_geom(ITAxes);
		  free(sset_dist);
	          break;

	 default: /*General case*/
/*	          sset_dist = init_array(nspher_set*(nspher_set+1)/2);
		  calc_distance(sset_geom,sset_dist,nspher_set);
		  a1 = 0;
		  a2 = (inv_related(sset_geom[a1],sset_geom[a1+1])) ? a1+2 : a1+1;
		  r = sset_dist[ioff[a2]+i1];
		  num_a3 = 0;
		  for(i=a2+1;i<nspher_set;i++)
		    if (fabs(r - sset_dist[ioff[i]+a1]) < ZERO) {
		      a3[num_a3] = i;
		      num_a3++;
		    }
		  switch (num_a3 + 1) {
		    case 1:
		    case 2: for(i=a2+1;i<nspher_set;i++)
			      if () {
			        a3[num_a3] = i;
				num_a3++;
			      }
			      */
	          break;
	}
	free_matrix(sset_geom,num_atoms);
      }
    }
    else if (num_atoms == 1) { /* Atom */
      fprintf(outfile,"    It is a spherical top.\n");
      rotor = atom;
    }
    else if (num_atoms <= 0) { /* ??? */
      punt("Fewer than 1 atom");
    }

    free(v1);
    free(v2);
    free(v3);
    free(IM);
    free_block(IT);
    free_block(ITAxes);
}

