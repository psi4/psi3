/*! \file
    \ingroup OPTKING
    \brief SIMPLES.CC member functions for simples class
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"

namespace psi { namespace optking {

extern void print_mat2(double **, int, int, FILE *);

void simples_class :: locate_id(int id, Intco_type *itype, int *sub_index, int *sub_index2) const {
  int I, i, found = 0;
  char error[100];

  for (i=0; i<stre.size(); ++i)
    if (stre[i].id == id) {
      ++found; *itype = STRE_TYPE; *sub_index = i;
    }

  for (i=0; i<bend.size(); ++i)
    if (bend[i].id == id) {
      ++found; *itype = BEND_TYPE; *sub_index = i;
    }

  for (i=0; i<tors.size(); ++i)
    if (tors[i].id == id) {
      ++found; *itype = TORS_TYPE; *sub_index = i;
    }

  for (i=0; i<out.size(); ++i)
    if (out[i].id == id) {
      ++found; *itype = OUT_TYPE; *sub_index = i;
    }

  for (i=0;i<linb.size();++i)
    if (linb[i].id == id) {
      ++found; *itype = LINB_TYPE; *sub_index = i;
    }

  for (i=0;i<frag.size();++i) {
    for (I=0; I<6; ++I) {
      if (frag[i].coord_on[I]) {
        if((frag[i].id + I) == id) {
          ++found; *itype = FRAG_TYPE; *sub_index = i; *sub_index2 = I;
        }
      }
    }
  }

  if (found == 0) {
    sprintf(error,"Error: Simple internal id %d not found.\n", id);
    throw(error);
  }
  else if (found > 1) {
    sprintf(error,"Error: Simple internal id %d appears more than once.\n", id);
    throw(error);
  }
  // fprintf(outfile,"returning intco_type: %d sub_index %d\n",*intco_type, *sub_index);
  return;
}


double ** simples_class::bond_connectivity_matrix(int natoms) const {
  int i, a, b;
  double **B = block_matrix(natoms,natoms);
  for (i=0; i<stre.size(); ++i) {
    a = stre[i].A;
    b = stre[i].B;
    B[b][a] = B[a][b] = 1;
  }
  return B;
}


/** Converts the absolute optking index into the user-assigned id number ***/
// skips over interfragment coordinates turned off
int simples_class :: index_to_id(int index) const {
  int i, I, count = -1, id = -1;

  for (i=0; i<stre.size(); ++i)
    if (++count == index) id = stre[i].id;
  for (i=0; i<bend.size(); ++i)
    if (++count == index) id = bend[i].id;
  for (i=0; i<tors.size(); ++i)
    if (++count == index) id = tors[i].id;
  for (i=0; i<out.size(); ++i)
    if (++count == index) id = out[i].id;
  for (i=0; i<linb.size(); ++i)
    if (++count == index) id = linb[i].id;

  /* we assume an id number of the fragment set plus the index 0-5 of the coordinate in the set */
  for (i=0; i<frag.size(); ++i) {
    for (I=0;I<6;++I)
      if (frag[i].get_coord_on(I))
        if (++count == index)
          id = frag[i].id + I;
  }

  // fprintf(outfile,"index_to_id(%d): returning id = %d\n",index, id);
  return id;
}

// converts id number to total optking index
// only interfragment coodinates turned on get counted in optking index
int simples_class :: id_to_index(int id) const {
  int i, I, count = 0, index = -1, found = 0;
  char error[100];

  for (i=0; i<stre.size(); ++i,++count)
    if (id == stre[i].id) {index = count; ++found;}

  for (i=0; i<bend.size(); ++i,++count)
    if (id == bend[i].id) {index = count; ++found;}

  for (i=0; i<tors.size(); ++i,++count)
    if (id == tors[i].id) {index = count; ++found;}

  for (i=0; i<out.size(); ++i,++count)
    if (id == out[i].id) {index = count; ++found;}

  for (i=0; i<linb.size(); ++i,++count)
    if (id == linb[i].id) {index = count; ++found;}

  /* we assume an id number of the fragment set plus the index 0-5 of the coordinate in the set */
  for (i=0; i<frag.size();++i) {
    for (I=0; I<6; ++I) {
      if (frag[i].get_coord_on(I)) {
        if (id == (frag[i].id + I)) {index = count; ++found;}
        ++count;
      }
    }
  }

  if (found == 0) {
    sprintf(error,"Error: Simple internal id %d not found.\n", id);
    throw(error);
  }
  else if (found > 1) {
    sprintf(error,"Error: Simple internal id %d appears more than once.\n", id);
    throw(error);
  }
  //  fprintf(outfile,"id_to_index(%d): returning index = %d\n",id, index);
  return index;
}

int simples_class::get_id_from_atoms_stre(int a, int b) const {
  int i;
  char error[100];

  stre_class s1(0, a, b);
  for (i=0; i<stre.size(); ++i) {
    if (s1 == stre[i])
      return (stre[i].id);
  }
  sprintf(error,"Could not find simple stretch for atoms %d %d", a+1, b+1);
  throw(error);
}

int simples_class::get_id_from_atoms_bend(int a, int b, int c) const {
  int i;
  char error[100];

  bend_class b1(0, a, b, c);
  for (i=0; i<bend.size(); ++i) {
    if (b1 == bend[i]) 
      return (bend[i].id);
  }
  sprintf(error,"Could not find simple bend for atoms %d %d %d", a+1, b+1, c+1);
  throw(error);
}

// multiplies sign by -1 if matching torsion is D-C-B-A
int simples_class::get_id_from_atoms_tors(int a, int b, int c, int d) const {
  int i;
  char error[100];

  tors_class t1(0, a, b, c, d);
  for (i=0; i<tors.size(); ++i) {
    if (t1 == tors[i])
      return (tors[i].id);
  }
  sprintf(error,"Could not find simple torsion for atoms %d %d %d %d", a+1, b+1, c+1, d+1);
  throw(error);
}

// multiplies sign by -1 if matching out-of-plane coordinate has C and D reversed
int simples_class::get_id_from_atoms_out(int a, int b, int c, int d, int *sign) const {
  int i;
  char error[100];

  for (i=0; i<out.size(); ++i) {
    if ( (out[i].A == a) && (out[i].B == b) && (out[i].C == c) && (out[i].D == d)) {
      return (out[i].id);
    }
    else if ( (out[i].A == a) && (out[i].B == b) && (out[i].D == c) && (out[i].C == d)) {
      *sign *= -1;
      return (out[i].id);
    }
  }
  sprintf(error,"Could not find simple out of plane for atoms %d %d %d %d", a+1, b+1, c+1, d+1);
  throw(error);
}

int simples_class::get_id_from_atoms_linb(int a, int b, int c, int linval) const {
  int i;
  char error[100];

  linb_class b1(0, a, b, c, linval);
  for (i=0; i<linb.size(); ++i) {
    if (b1 == linb[i])
      return (linb[i].id);
  }
  sprintf(error,"Could not find simple linb for atoms %d %d %d", a+1, b+1, c+1);
  throw(error);
}

int simples_class::get_id_from_atoms_frag(int a_natom, int b_natom, int *a_atom, int *b_atom) const {
  int i, a, b, match=0;
  char error[100];

  while ((match==0) && (i<frag.size())) {
    match = 1;
    if (a_natom != frag[i].A_natom) match = 0;
    if (b_natom != frag[i].B_natom) match = 0;
    if (match) {
      for (a=0; a<a_natom; ++a) {
        if (a_atom[a] != frag[i].A_atom[a])
          match=0;
      }
      for (b=0; b<b_natom; ++b) {
        if (b_atom[b] != frag[i].B_atom[b])
          match=0;
      }
    }
    if (match)
      return frag[i].id;
    ++i;
    // lets return id number of first element of set - R(AB) for now
  }

  sprintf(error,"Could not find simple fragment with natoms %d %d.\n", a_natom, b_natom);
  throw(error);
}


/** constructs simples_class from ip parsing or geometry **/
simples_class :: simples_class(cartesians& carts, int user_intcos)
           : stre(), bend(), tors(), out(), linb(), frag() {

    int id, i,I,j,k,a,b,c,d,e,lin,dim,A_P,B_P,n;
    int na, nb;
    double *v, dot;
    char error[100];

    stre_class *s1;
    bend_class *b1;
    tors_class *t1;
    out_class *o1;
    linb_class *lb1;
    frag_class *f1;

    if (user_intcos) { // Read in simple internal coordinates from intco.dat
    try {

      if (ip_exist("STRE",0)) {
        n=0;
        ip_count("STRE",&n,0);

        for(i=0;i<n;++i) {
          ip_count("STRE",&j,1,i);
          if (j != 3) {
            sprintf(error,"Stretch %d is of wrong dimension.\n",i+1);
            throw(error);
          }
          ip_data("STRE", "%d", &id, 2, i, 0);
          ip_data("STRE", "%d",  &a, 2, i, 1);
          ip_data("STRE", "%d",  &b, 2, i, 2);

          s1 = new stre_class(id, a-1, b-1);
          stre.push_back(*s1);
        }
      }

      if (ip_exist("BEND",0)) {
        n=0;
        ip_count("BEND",&n,0);

        for(i=0;i<n;++i) {
          ip_count("BEND",&j,1,i);
          if (j != 4) {
            sprintf(error,"Bend %d is of wrong dimension.\n",i+1);
            throw(error);
          }
          ip_data("BEND", "%d", &id, 2, i, 0);
          ip_data("BEND", "%d",  &a, 2, i, 1);
          ip_data("BEND", "%d",  &b, 2, i, 2);
          ip_data("BEND", "%d",  &c, 2, i, 3);

          b1 = new bend_class(id, a-1, b-1, c-1);
          bend.push_back(*b1);
        }
      }

      if (ip_exist("LIN1",0)) {
        n=0;
        ip_count("LIN1",&n,0);
        for(i=0;i<n;++i) {
          ip_count("LIN1",&j,1,i);
          if (j != 4) {
            sprintf(error,"Linear bend %d is of wrong dimension.\n",i+1);
            throw(error);
          }
          ip_data("LIN1", "%d", &id, 2, i, 0);
          ip_data("LIN1", "%d",  &a, 2, i, 1);
          ip_data("LIN1", "%d",  &b, 2, i, 2);
          ip_data("LIN1", "%d",  &c, 2, i, 3);
          lb1 = new linb_class(id, a-1, b-1, c-1, 1);
          linb.push_back(*lb1);
        }
      }

      if (ip_exist("LIN2",0)) {
        n=0;
        ip_count("LIN2",&n,0);
        for(i=0;i<n;++i) {
          ip_count("LIN2",&j,1,i);
          if (j != 4) {
            sprintf(error,"Linear bend %d is of wrong dimension.\n",i+1);
            throw(error);
          }
          ip_data("LIN2", "%d", &id, 2, i, 0);
          ip_data("LIN2", "%d",  &a, 2, i, 1);
          ip_data("LIN2", "%d",  &b, 2, i, 2);
          ip_data("LIN2", "%d",  &c, 2, i, 3);
          lb1 = new linb_class(id, a-1, b-1, c-1, 2);
          linb.push_back(*lb1);
        }
      }

      if (ip_exist("TORS",0)) {
        n = 0;
        ip_count("TORS",&n,0);

        for(i=0;i<n;++i) {
          ip_count("TORS",&j,1,i);
          if (j != 5) {
            sprintf(error,"Torsion %d is of wrong dimension.\n",i+1);
            throw(error);
          }
          ip_data("TORS", "%d", &id, 2, i, 0);
          ip_data("TORS", "%d",  &a, 2, i, 1);
          ip_data("TORS", "%d",  &b, 2, i, 2);
          ip_data("TORS", "%d",  &c, 2, i, 3);
          ip_data("TORS", "%d",  &d, 2, i, 4);

          t1 = new tors_class(id, a-1, b-1, c-1, d-1);
          tors.push_back(*t1);
        }
      }

      if (ip_exist("OUT",0)) {
        n = 0;
        ip_count("OUT",&n,0);

        for(i=0; i<n; ++i) {
          ip_count("OUT",&j,1,i);
          if (j != 5) {
            sprintf(error,"Out-of-plane %d is of wrong dimension.\n",i+1);
            throw(error);
          }
          ip_data("OUT", "%d", &id, 2, i, 0);
          ip_data("OUT", "%d",  &a, 2, i, 1);
          ip_data("OUT", "%d",  &b, 2, i, 2);
          ip_data("OUT", "%d",  &c, 2, i, 3);
          ip_data("OUT", "%d",  &d, 2, i, 4);

          o1 = new out_class(id, a-1, b-1, c-1, d-1);
          out.push_back(*o1);
        }
      }

      if (ip_exist("FRAG",0)) {
        n = 0;
        ip_count("FRAG",&n,0);
        for(i=0; i<n; ++i) {
          // id, coord_on, atoms_A, atoms_B, ref A weights, ref B weights
          ip_count("FRAG",&j,1,i);
          if (j != 6) {
            sprintf(error,"Fragment %d should have 6 parts to its definiction",i+1);
            throw(error);
          }

          ip_data("FRAG","%d", &id, 2, i, 0);
          ip_count("FRAG",     &na, 2, i, 2);
          ip_count("FRAG",     &nb, 2, i, 3);

          // make sure 3 reference points are specified for each fragment (for now)
          ip_count("FRAG", &A_P, 2, i, 4);
          ip_count("FRAG", &B_P, 2, i, 5);
          if ( (A_P != 3) || (B_P != 3)) {
            sprintf(error,"Fragment %d should give weights for 3 reference points for each fragment", i+1);
            throw(error);
          }

          f1 = new frag_class(id, na, nb, A_P, B_P);
          frag.push_back(*f1);

          // read what coordinates are on
          ip_count("FRAG",&a,2,i,1);
          if (a != 6)
            throw("Coordinate-on (second) part of Fragment %d must have dimension 6.\n",i+1);
          for (dim=0, I=0; I<6;++I) {
            ip_boolean("FRAG",&a,3,i,1,I);
            frag[i].set_coord_on(I,a); // set coordinate on or off
            if (a) ++dim;
          }
          // if only interfragment stretch, then use only 1 reference atom
          ip_boolean("FRAG",&a,3,i,1,0);
          if (a && (dim == 1)) {
            A_P = B_P = 1;
            frag[i].set_A_P(A_P);
            frag[i].set_B_P(B_P);
          }

          // read what atoms are in each fragment
          for (j=0;j<na;++j) {
            ip_data("FRAG","%d",&a,3,i,2,j);
            frag[i].set_A_atom(j,a-1);  
          }
          for (j=0;j<nb;++j) {
            ip_data("FRAG","%d",&b,3,i,3,j);
            frag[i].set_B_atom(j,b-1); 
          }

          // read in weights for fragment A - make the absolute sum of the coefficients 1
          v = new double[na];
          for (j=0;j<A_P;++j) {
            ip_count("FRAG",&a,3,i,4,j);
            if (a != na) throw("Fragment A has weights of wrong dimension\n");
            for (a=0;a<na;++a)
              ip_data("FRAG","%lf",&(v[a]),4,i,4,j,a);
            dot = 0.0;
            for (a=0;a<na;++a)
              dot += fabs(v[a]);
            for (a=0;a<na;++a)
              v[a] /= dot;
            for (a=0;a<na;++a)
              frag[i].set_A_weight(j,a,v[a]);
          }
          delete [] v;

          // read in weights for fragment B - make the absolute sum of the coefficients 1
          v = new double[nb];
          for (j=0;j<B_P;++j) {
            ip_count("FRAG",&b,3,i,5,j);
            if (b != nb) throw("Fragment A has weights of wrong dimension\n");
            for (b=0;b<nb;++b)
              ip_data("FRAG","%lf",&(v[b]),4,i,5,j,b);
            dot = 0.0;
            for (b=0;b<nb;++b)
              dot += fabs(v[b]);
            for (b=0;b<nb;++b)
              v[b] /= dot;
            for (b=0;b<nb;++b)
              frag[i].set_B_weight(j,b,v[b]);
          }
          delete [] v;
        }
      }
    }
    catch(char const *str) {
      fprintf(outfile,"Trouble reading simple internal coordinates\n");
      fprintf(stderr,"Trouble reading simple internal coordinates\n");
      punt(str);
    }
    }
    else { // Generate simple internal coordinates from a cartesian geometry
    try {
      int count, Z1, Z2, natom;
      double **atom_dist, **coord_2d, tval;
      int **bonds, id_count=0;
      int num_bonds, num_nobonds, n_connected, n_connected_old;

      natom = carts.get_natom();
      coord_2d = carts.get_coord_2d();

      fprintf(outfile,"\nGenerating simple internals\n"); fflush(outfile);

      /* Compute atomic distance matrix */
      atom_dist = block_matrix(natom,natom);
      for (i=0; i<natom; ++i)
        for (j=0; j<natom; ++j)
          atom_dist[i][j] = sqrt(SQR(coord_2d[i][0] - coord_2d[j][0])+
                                 SQR(coord_2d[i][1] - coord_2d[j][1])+
                                 SQR(coord_2d[i][2] - coord_2d[j][2]));

      /* Determine bond connectivity matrix using distance criteria */
      bonds = init_int_matrix(natom, natom);
      for (i=0; i<natom; ++i) {
        Z1 = (int) carts.get_Z(i);
        if (Z1 == 0) continue;
        for (j=0; j<i; ++j) {
          Z2 = (int) carts.get_Z(j);
          if (Z2 == 0) continue;
          if (cov_radii[Z2] != 0.0) {
            tval = (cov_radii[Z1] + cov_radii[Z2])/_bohr2angstroms;// to au
            if (atom_dist[i][j] < (optinfo.scale_connectivity * tval))
              bonds[i][j] = bonds[j][i] = 1;
          }
        }
      }

      /* Add user-defined bonds and nobonds */
      if (!ip_exist("BONDS",0))
        num_bonds = 0;
      else
        ip_count("BONDS",&num_bonds,0,0);
      if (optinfo.print_simples)
        fprintf(outfile,"\nNumber of user-specified bonds: %3d\n",num_bonds);
      for(i=0;i<num_bonds;++i) {
        ip_data("BONDS","%d",&a,2,i,0);
        ip_data("BONDS","%d",&b,2,i,1);
        fprintf(outfile,"Bond %3d: %3d %3d\n",i+1,a,b);
        a -= 1;
        b -= 1;
        bonds[a][b] = 1;
        bonds[b][a] = 1;
      }

      if (!ip_exist("NOBONDS",0))
        num_nobonds = 0;
      else
        ip_count("NOBONDS",&num_nobonds,0,0);
      if (optinfo.print_simples)
        fprintf(outfile,"\nNumber of user-specified \"nobonds\": %3d\n",num_nobonds);
      for(i=0;i<num_nobonds;++i) {
        ip_data("NOBONDS","%d",&a,2,i,0);
        ip_data("NOBONDS","%d",&b,2,i,1);
        fprintf(outfile,"Nobond %3d: %3d %3d\n",i+1,a,b);
        a -= 1;
        b -= 1;
        bonds[a][b] = 0;
        bonds[b][a] = 0;
      }

      // if we are taking a geometry step or displacing along internal coordinates, then
      // check to make sure that all atoms are bonded
      if (optinfo.mode == MODE_OPT_STEP || optinfo.mode == MODE_DISP_IRREP ||
          optinfo.mode == MODE_DISP_LOAD) {
        int *atoms_connected_to_first = new int[natom];
        for (i=0;i<natom;++i)
          atoms_connected_to_first[i] = 0;
        n_connected = 1;
        atoms_connected_to_first[0] = 1;
      
        do {
          n_connected_old = n_connected;
          for (i=0;i<natom;++i) {
            if (atoms_connected_to_first[i]) {
              for (j=0;j<natom;++j) {
                if (!atoms_connected_to_first[j] && bonds[i][j]) { // a new connection
                  ++n_connected;
                  atoms_connected_to_first[j] = 1;
                }
              }
            }
          }
        }
        while (n_connected > n_connected_old);
  
        n_connected = 0;
        for (i=0;i<natom;++i)
          if (atoms_connected_to_first[i]) ++n_connected;
        
        if (n_connected != natom) {
          fprintf(outfile,"\n\t** ERROR - Not all atoms are connected by bonds.\n");
          fprintf(outfile,"\tYou will need to specify internal coordinates manually.\n");
          fprintf(outfile,"\n\tAtoms connected to the 1st atom:");
          for (i=1;i<natom;++i)
            if (atoms_connected_to_first[i]) fprintf(outfile," %d",i+1);
          fprintf(outfile,"\n");
          abort();
        }
        delete [] atoms_connected_to_first;
      } // end - check to make sure that all atoms are bonded

      /* determine stretch internal coordinates */
      id_count = 0; // index the id numbers of simples
      for (i=0; i<natom; ++i)
        for (j=i+1; j<natom; ++j)
          if (bonds[i][j]) {  // smaller atom first convention
            s1 = new stre_class(++id_count, i, j);
            stre.push_back(*s1);
          }

      if (optinfo.print_simples) {
        fprintf(outfile,"\n+++ Bond Connectivity +++\n");
        for (i=0; i<natom; ++i) {
          fprintf(outfile,"%d:",i+1);
          for (j=0; j<natom; ++j)
            if (bonds[i][j]) fprintf(outfile," %d",j+1);
          fprintf(outfile,"\n");
        }
      }

      /* Determine bond angles */
      for(i=0;i<natom;++i)
        for(j=0;j<natom;++j)
          for(k=i+1; k<natom; ++k)
            if (bonds[i][j] && bonds[j][k]) {
              /* check if angle is linear */
              tval = atom_dist[i][k] - atom_dist[i][j] - atom_dist[j][k];
              if ( fabs(tval) > NONLINEAR_DIST ) {
                b1 = new bend_class(++id_count, i, j, k);
                bend.push_back(*b1);
              }
            }

      /* Determine torsions */
      for(a=0;a<natom;++a)
        for(b=0;b<natom;++b)
          if (a!=b) {
            for(c=0;c<natom;++c)
              if ( (c!=a) && (c!=b) ) {
                for(d=a+1; d<natom; ++d)
                  if ((d!=b) && (d!=c) && bonds[a][b] && bonds[b][c] && bonds[c][d]) {
                    /* don't include if there are linear angles a-b-c or b-c-d */
                    tval = atom_dist[a][c] - atom_dist[a][b] - atom_dist[b][c];
                    lin = 0;
                    if ( fabs(tval) < NONLINEAR_DIST ) lin = 1;
                    tval = atom_dist[b][d] - atom_dist[b][c] - atom_dist[c][d];
                    if ( fabs(tval) < NONLINEAR_DIST ) lin = 1;
                    if (!lin) {
                      t1 = new tors_class(++id_count, a, b, c, d);
                      tors.push_back(*t1);
                    }
                  }
              }
          }

      // add bonus torsions aa-a-b-c-cc if a-b-c is linear, etc.
      int aa, bb, cc, ia, ib, ic, ba, bc, lin_aa, lin_cc, found_aa, found_cc;
      int *taa, *tbb, *tcc, *tdd, skip;
      taa = new int[300]; tbb = new int[300]; tcc = new int[300]; tdd = new int[300];

      for(ia=0;ia<natom;++ia)
        for(ib=0;ib<natom;++ib) if ( (ia!=ib) && bonds[ia][ib] )
          for(ic=0;ic<natom;++ic) if ( (ic!=ib) && bonds[ib][ic] ) {

            a = ia; b = ib; c = ic;
            tval = atom_dist[a][c] - atom_dist[a][b] - atom_dist[b][c];
            if ( fabs(tval) < NONLINEAR_DIST ) {
            // a-b-c is linear

          // if ib is bonded to any other atoms, then ignore this linear segment
          for (skip=0,bb=0;bb<natom;++bb)
            if (bb!=ia && bb!=ic && bonds[bb][ib]) skip = 1;
          if (skip) continue; // continue to next ic

              fprintf(outfile," %d %d %d is colinear\n",a,b,c);
 
              // find aa bonded to a such that aa-a-b is not linear too
              found_aa = 0;
              lin_aa = a;
              a = b; // atom adjacent to a
              while (!found_aa) {
              ba = a;
              a = lin_aa;
              for (aa=0;aa<natom;++aa) if (aa!=ba && bonds[aa][a]) {
                c = ic;
                tval = atom_dist[aa][ba] - atom_dist[aa][a] - atom_dist[a][ba];
                if ( fabs(tval) > NONLINEAR_DIST ) { 
                // aa-a-b is not linear
                  found_aa = 1;
                  // fprintf(outfile,"found aa to a,b,c to %d, %d %d %d\n", aa,a,ba,c);

                  // find cc bonded to c such that b-c-cc is not linear or b-c-d-cc,etc.
                  found_cc = 0;
                  lin_cc = c;
                  c = b; // atom adjacent to c
                  while (!found_cc) {
                  bc = c;
                  c = lin_cc;
                  for (cc=aa+1;cc<natom;++cc) if (cc!=bc && bonds[cc][c]) {
                    tval = atom_dist[cc][bc] - atom_dist[cc][c] - atom_dist[bc][c];
                    if ( fabs(tval) > NONLINEAR_DIST ) { 
                    // b-c-cc is not linear
                      // fprintf(outfile,"found cc to a,b,c to %d %d %d, %d\n", a,bc,c,cc);
                      found_cc = 1;
                      t1 = new tors_class(0, aa, a, c, cc);
                      if (is_unique(*t1)) {
                        delete t1;
                        t1 = new tors_class(++id_count, aa, a, c, cc);
                        tors.push_back(*t1);
                      }
                      else
                        delete t1;
                    }
                    else lin_cc = cc;
                  } // end loop over cc
                  if ( !found_cc && (c == lin_cc) ) break;
                  } // end while search for cc bonded to c
                } // end loop over aa-a-b is not linear
                else lin_aa = aa;
              } // end loop over aa
              if ( !found_aa && (a == lin_aa) ) break;
              }
            } // end loop over linear a-b-c
          } // end loop over ia, ib, ic

      delete taa, tbb, tcc, tdd;

      // add bonus torsions for cases like BH3 and H2CO
      // where natom = 3 and total number of bonds is 3
      if (natom == 4) {
        int nbonds=0;
        for (a=0;a<natom;++a)
          for (b=0;b<a;++b) {
            if (bonds[a][b])
              nbonds++;
          }
        if (nbonds == 3) {
          // make sure whole molecule is not linear
          if (bend.size() > 2) {
            // we'll just try 3 torsions for now - this works with NH3 with different
            // atom numberings -- these torsions may have to be changed to
            // use the central atom in a specific way someday
            t1 = new tors_class(++id_count, 0, 1, 2, 3, 0.0);
            tors.push_back(*t1);

            t1 = new tors_class(++id_count, 0, 2, 3, 1, 0.0);
            tors.push_back(*t1);

            t1 = new tors_class(++id_count, 0, 3, 1, 2, 0.0);
            tors.push_back(*t1);
          }
        }
      }

      free_block(coord_2d);
      free_block(atom_dist);

      ffile(&fp_intco, "intco.dat",0);
      print(fp_intco,0); 
      fclose(fp_intco);
      
      free_int_matrix(bonds);
    }
    catch (const char *str) {
      fprintf(outfile,"Trouble generating simple internal coordinates\n");
      fprintf(stderr,"Trouble generating simple internal coordinates\n");
      punt(str);
    }
    } // end generate simples

    if (get_num() == 0) {
      fprintf(outfile,"You may need to increase the value of the scale_connectivity keyword.");
      sprintf(error, "Error: No simple internals were generated and read.");
      throw(error);
    }
    return;
  }


void simples_class :: print(FILE *fp_out, bool print_vals, bool print_frag_weights) const {
  int i, I,a,b;

  if (!print_vals) fprintf(fp_out, "intco: (\n");

  if (stre.size() > 0) {
    if (!print_vals) fprintf(fp_out,"  stre = (\n");
    else fprintf(fp_out,"Stretches\n");
    for (i=0; i < stre.size(); ++i)
      stre[i].print(fp_out, print_vals);
    if (!print_vals) fprintf(fp_out,"  )\n");
  }

  if (bend.size() > 0) {
    if (!print_vals) fprintf(fp_out,"  bend = (\n");
    else fprintf(fp_out, "Bends\n");
    for (i=0; i < bend.size(); ++i)
      bend[i].print(fp_out, print_vals);
    if (!print_vals) fprintf(fp_out,"  )\n");
  }

  if (tors.size() > 0) {
    if (!print_vals) fprintf(fp_out,"  tors = (\n");
    else fprintf(fp_out,"Torsions\n");
    for (i=0; i < tors.size(); ++i)
      tors[i].print(fp_out, print_vals);
    if (!print_vals) fprintf(fp_out,"  )\n");
  }

  if (out.size() > 0) {
    if (!print_vals) fprintf(fp_out,"  out = (\n");
    else fprintf (fp_out,"Out-of-planes\n");
    for (i=0; i < out.size(); ++i)
      out[i].print(fp_out, print_vals);
    if (!print_vals) fprintf(fp_out,"  )\n");
  }

  if (linb.size() > 0) {
    if (!print_vals) {
      fprintf(fp_out,"  lin1 = (\n");
      for (i=0; i < linb.size(); ++i) {
        if (linb[i].get_linval() == 1)
          linb[i].print(fp_out, print_vals);
      }
      fprintf(fp_out,"  )\n");
      fprintf(fp_out,"  lin2 = (\n");
      for (i=0; i < linb.size(); ++i) {
        if (linb[i].get_linval() == 2)
          linb[i].print(fp_out, print_vals);
      }
      fprintf(fp_out,"  )\n");
    }
    else {
      fprintf(fp_out, "Linear Bends\n");
      for (i=0; i < linb.size(); ++i)
        linb[i].print(fp_out, print_vals);
    }
  }

  if (frag.size() > 0) {
    if (!print_vals) fprintf(fp_out,"  fragment = (\n");
    else fprintf(fp_out, "Fragments\n");
    for (i=0; i < frag.size(); ++i)
      frag[i].print(fp_out, print_vals, print_frag_weights);
    if (!print_vals) fprintf(fp_out,"  )\n");
  }

  if (!print_vals) fprintf(fp_out, ")\n");
  return;
}

}} /* namespace psi::optking */

