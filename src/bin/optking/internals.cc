/*** INTERNALS.CC member functions for internals class ***/ 

#if HAVE_CMATH
# include <cmath>
#else
# include <math.h>
#endif

extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <physconst.h>
#include "bond_lengths.h"
}

#define EXTERN
#include "opt.h"
#undef EXTERN
#include "cartesians.h"
#include "internals.h"

int is_unique(int a, int b, int c, int d, int *aa, int *bb, int *cc, int *dd);

/*** INTERNALS::INDEX_TO_ID
 * Converts the absolute optking index into the user-assigned id number ***/

int internals :: index_to_id(int index) {
  int i, count = -1, id = -1;
  for (i=0;i<stre.get_num();++i)
    if (++count == index) {id = stre.get_id(i); break;}
  for (i=0;i<bend.get_num();++i)
    if (++count == index) {id = bend.get_id(i); break;}
  for (i=0;i<tors.get_num();++i)
    if (++count == index) {id = tors.get_id(i); break;}
  for (i=0;i<out.get_num();++i)
    if (++count == index) {id = out.get_id(i); break;}
  for (i=0;i<lin_bend.get_num();++i)
    if (++count == index) {id = lin_bend.get_id(i); break;}
  //   fprintf(outfile,"index_to_id(%d): returning id = %d\n",index, id);
  return id;
}



/*** INTERNALS::ID_TO_INDEX
 * Converts the user-assigned id number into an absolute optking index ***/ 

int internals :: id_to_index(int id) {
  int i, count = 0, index = -1, found = 0;
  for (i=0;i<stre.get_num();++i,++count)
    if (id == stre.get_id(i)) {index = count; ++found;}
  for (i=0;i<bend.get_num();++i,++count)
    if (id == bend.get_id(i)) {index = count; ++found;}
  for (i=0;i<tors.get_num();++i,++count)
    if (id == tors.get_id(i)) {index = count; ++found;}
  for (i=0;i<out.get_num();++i,++count)
    if (id == out.get_id(i)) {index = count; ++found;}
  for (i=0;i<lin_bend.get_num();++i,++count)
    if (id == lin_bend.get_id(i)) {index = count; ++found;}
  if (found == 0) {
    fprintf(outfile,"Error: Simple internal id not found.\n");
    exit(2);
  }
  else if (found > 1) {
    fprintf(outfile,"Error: Simple internal id appears more than once.\n");
    exit(2);
  }
  // fprintf(outfile,"id_to_index(%d): returning index = %d\n",id, index);
  return index;
}



/*** INTERNALS::LOCATE_ID
 * Given the user-assigned id number, returns the type of internal
 and the subindex within that type ***/

void internals :: locate_id(int id, int *intco_type, int *sub_index) {
  int i, found = 0;
  // fprintf(outfile,"locating id %d\n",id);
  for (i=0;i<stre.get_num();++i)
    if (stre.get_id(i) == id) {
      ++found; *intco_type = STRE_TYPE; *sub_index = i;
    }
  for (i=0;i<bend.get_num();++i)
    if (bend.get_id(i) == id) {
      ++found; *intco_type = BEND_TYPE; *sub_index = i;
    }
  for (i=0;i<tors.get_num();++i)
    if (tors.get_id(i) == id) {
      ++found; *intco_type = TORS_TYPE; *sub_index = i;
    }
  for (i=0;i<out.get_num();++i)
    if (out.get_id(i) == id) {
      ++found; *intco_type = OUT_TYPE; *sub_index = i;
    }
  for (i=0;i<lin_bend.get_num();++i)
    if (lin_bend.get_id(i) == id) {
      ++found; *intco_type = LIN_BEND_TYPE; *sub_index = i;
    }
  if (found == 0) {
    fprintf(outfile,"Error: Simple internal id not found.\n");
    exit(2);
  }
  else if (found > 1) {
    fprintf(outfile,"Error: Simple internal id appears more than once.\n");
    exit(2);
  }
  // fprintf(outfile,"returning intco_type: %d sub_index %d\n",*intco_type, *sub_index);
  return;
}

internals :: ~internals() {
  // fprintf(stdout,"destructing internals\n");
  // Destructors for stre, bend, tors, out are called automatically
}

// minimal constructor
internals :: internals(int *size_arr) : stre(size_arr[0]), bend(size_arr[1]),
  tors(size_arr[2]), out(size_arr[3]), lin_bend(size_arr[4]) {
}

/*** INTERNALS::INTERNALS class constructor for internals ***/ 

internals :: internals(cartesians& carts, int user_intcos)
           : stre(), bend(), tors(), out(), lin_bend() {

    int i,j,k,a,b,c,d,e,cnt,lin, type_count=0;
    int cnt_stre, cnt_bend, cnt_tors, cnt_out, cnt_lin_bend;
    int cnt_lin1_bend, cnt_lin2_bend;
    cnt_stre = cnt_bend = cnt_tors = cnt_out = cnt_lin_bend = 0;
    cnt_lin1_bend = cnt_lin2_bend = 0;

    for (cnt=1; cnt>=0; --cnt) { /* 1st time through we are just counting */

      if (!cnt) { /* second time through allocate memory and add internals */
        stre.set_num(cnt_stre);
        bend.set_num(cnt_bend);
        tors.set_num(cnt_tors);
        out.set_num(cnt_out);
        lin_bend.set_num(cnt_lin_bend);

        stre.allocate(cnt_stre);
        bend.allocate(cnt_bend);
        tors.allocate(cnt_tors);
        out.allocate(cnt_out);
        lin_bend.allocate(cnt_lin_bend);
      }

    if (user_intcos) {
      /* Read in simple internal coordinates from intco.dat */

      cnt_stre = 0;
      if (ip_exist("STRE",0)) {
        i=0;
        ip_count("STRE",&i,0);
        cnt_stre = i;

        for(i=0;i<cnt_stre;++i) {
          ip_count("STRE",&j,1,i);
          if (j != 3) {
            fprintf(outfile,"Stretch %d is of wrong dimension.\n",i+1);
            exit(2);
          }
          ip_data("STRE","%d",&(a),2,i,0);
          ip_data("STRE","%d",&(b),2,i,1);
          ip_data("STRE","%d",&(c),2,i,2);
          if (!cnt) {
            stre.set_id(i,a);
            swap(&b,&c);
            stre.set_A(i,b-1);
            stre.set_B(i,c-1);
            stre.set_val(i,0.0);
          }
        }
      }

      cnt_bend = 0;
      if (ip_exist("BEND",0)) {
        i = 0;
        ip_count("BEND",&i,0);
        cnt_bend = i;
        for(i=0;i<cnt_bend;++i) {
          ip_count("BEND",&j,1,i);
          if (j != 4) {
            fprintf(outfile,"Bend %d is of wrong dimension.\n",i+1);
            exit(2);
          }
          ip_data("BEND","%d",&(a),2,i,0);  ip_data("BEND","%d",&(b),2,i,1);
          ip_data("BEND","%d",&(c),2,i,2);  ip_data("BEND","%d",&(d),2,i,3);
          if (!cnt) {
            bend.set_id(i,a);
            swap(&b,&d);
            bend.set_A(i,b-1);
            bend.set_B(i,c-1);
            bend.set_C(i,d-1);
            bend.set_val(i,0.0);
          }
        }
      }

      // lin1 bends
      type_count = 0;
      if (ip_exist("LIN1",0)) {
        cnt_lin1_bend = 0;
        ip_count("LIN1",&cnt_lin1_bend,0);
        for(i=0;i<cnt_lin1_bend;++i) {
          ip_count("LIN1",&j,1,i);
          if (j != 4) {
            fprintf(outfile,"Linear bend %d is of wrong dimension.\n",i+1);
            exit(2);
          }
          ip_data("LIN1","%d",&(a),2,i,0);  ip_data("LIN1","%d",&(b),2,i,1);
          ip_data("LIN1","%d",&(c),2,i,2);  ip_data("LIN1","%d",&(d),2,i,3);
          if (!cnt) {
            lin_bend.set_id(type_count,a);
            swap(&b,&d);
            lin_bend.set_A(type_count,b-1);
            lin_bend.set_B(type_count,c-1);
            lin_bend.set_C(type_count,d-1);
            lin_bend.set_linval(type_count,1);
            lin_bend.set_val(type_count++,0.0);
          }
        }
      }
      if (ip_exist("LIN2",0)) {
        cnt_lin2_bend = 0;
        ip_count("LIN2",&cnt_lin2_bend,0);
        for(i=0;i<cnt_lin2_bend;++i) {
          ip_count("LIN2",&j,1,i);
          if (j != 4) {
            fprintf(outfile,"Linear bend %d is of wrong dimension.\n",i+1);
            exit(2);
          }
          ip_data("LIN2","%d",&(a),2,i,0);  ip_data("LIN2","%d",&(b),2,i,1);
          ip_data("LIN2","%d",&(c),2,i,2);  ip_data("LIN2","%d",&(d),2,i,3);
          if (!cnt) {
            lin_bend.set_id(type_count,a);
            swap(&b,&d);
            lin_bend.set_A(type_count,b-1);
            lin_bend.set_B(type_count,c-1);
            lin_bend.set_C(type_count,d-1);
            lin_bend.set_linval(type_count,2);
            lin_bend.set_val(type_count++,0.0);
          }
        }
      }
      cnt_lin_bend = cnt_lin1_bend + cnt_lin2_bend;

      cnt_tors = 0;
      if (ip_exist("TORS",0)) {
        i = 0;
        ip_count("TORS",&i,0);
        cnt_tors = i;

        for(i=0;i<cnt_tors;++i) {
          ip_count("TORS",&j,1,i);
          if (j != 5) {
            fprintf(outfile,"Torsion %d is of wrong dimension.\n",i+1);
            exit(2);
          }
          ip_data("TORS","%d",&(a),2,i,0);
          ip_data("TORS","%d",&(b),2,i,1);
          ip_data("TORS","%d",&(c),2,i,2);
          ip_data("TORS","%d",&(d),2,i,3);
          ip_data("TORS","%d",&(e),2,i,4);
          if (!cnt) {
            tors.set_id(i,a);
            swap_tors(&b,&c,&d,&e);
            tors.set_A(i,b-1);
            tors.set_B(i,c-1);
            tors.set_C(i,d-1);
            tors.set_D(i,e-1);
            tors.set_val(i,0.0);
          }
        }
      }

      if (ip_exist("OUT",0)) {
        i = 0;
        ip_count("OUT",&i,0);
        cnt_out = i;

        for(i=0; i<cnt_out; ++i) {
          ip_count("OUT",&j,1,i);
          if (j != 5) {
            fprintf(outfile,"Out-of-plane %d is of wrong dimension.\n",i+1);
            exit(2);
          }
          ip_data("OUT","%d",&(a),2,i,0);
          ip_data("OUT","%d",&(b),2,i,1);
          ip_data("OUT","%d",&(c),2,i,2);
          ip_data("OUT","%d",&(d),2,i,3);
          ip_data("OUT","%d",&(e),2,i,4);
          if (!cnt) {
            out.set_id(i,a);
            out.set_A(i,b-1);
            out.set_B(i,c-1);
            out.set_C(i,d-1);
            out.set_D(i,e-1);
            out.set_val(i,0.0);
          }
        }
      }
      else if(!ip_exist("OUT",0))
        cnt_out = 0;
    }
    else {

      /* Generate simple internal coordinates from a cartesian geometry */
      int *ioff, count, Z1, Z2, natom;
      double **atom_dist, **coord_2d, tval;
      int **bonds, id_count=0;
      int num_bonds, num_nobonds;

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
        Z1 = (int) carts.get_fatomic_num(i);
        if (Z1 == 0) Z1 = 10; // let dummy atom be like C
        if ( radii[Z1] == 0) {
          fprintf(outfile,"WARNING! Optking does not know what bond lengths");
          fprintf(outfile,"to expect for atom %d.\n",i+1);
          fprintf(outfile,"You will have to specify this atom's connectivity in input.");
          continue;
        }
        for (j=0; j<i; ++j) {
          Z2 = (int) carts.get_fatomic_num(j);
          if (Z2 == 0) Z2 = 10; // let dummy atom be like C
          if (radii[Z2] != 0.0) {
            tval = (radii[Z1] + radii[Z2])/(100.0*_bohr2angstroms);// to au
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

      /* determine stretch internal coordinates */
      id_count = 0; // index the id numbers of internals
      type_count = 0;
      for (i=0; i<natom; ++i)
        for (j=i+1; j<natom; ++j)
          if (bonds[i][j]) {  // smaller atom first convention
            if (cnt)
              type_count++;
            else {
              stre.set_id(type_count,++id_count);
              stre.set_A(type_count,i);
              stre.set_B(type_count,j);
              stre.set_val(type_count++,0.0); // values determined later
            }
          }
      cnt_stre = type_count;

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
      type_count = 0;
      for(i=0;i<natom;++i)
        for(j=0;j<natom;++j)
          for(k=i+1; k<natom; ++k)
            if (bonds[i][j] && bonds[j][k]) {
              /* check if angle is linear */
                tval = atom_dist[i][k] - atom_dist[i][j] - atom_dist[j][k];
                lin = 0;
                if ( fabs(tval) < NONLINEAR_DIST ) lin = 1;
                if (!lin) {
                  if (cnt)
                    type_count++;
                  else {
                    bend.set_id(type_count,++id_count);
                    bend.set_A(type_count,i);
                    bend.set_B(type_count,j);
                    bend.set_C(type_count,k);
                    bend.set_val(type_count++,0.0);
                  }
                }
            }
      cnt_bend = type_count;

      /* Determine linear bond angles - code matches that above */
      type_count = 0;
      for(i=0;i<natom;++i)
        for(j=0;j<natom;++j)
          for(k=i+1; k<natom; ++k)
            if (bonds[i][j] && bonds[j][k]) {
              /* check if angle is linear */
                tval = atom_dist[i][k] - atom_dist[i][j] - atom_dist[j][k];
                lin = 0;
                if ( fabs(tval) < NONLINEAR_DIST ) {
                  if (cnt)
                    type_count += 2;
                  else {
                    lin_bend.set_id(type_count,++id_count);
                    lin_bend.set_A(type_count,i);
                    lin_bend.set_B(type_count,j);
                    lin_bend.set_C(type_count,k);
                    lin_bend.set_linval(type_count,1);
                    lin_bend.set_val(type_count++,0.0);

                    lin_bend.set_id(type_count,++id_count);
                    lin_bend.set_A(type_count,i);
                    lin_bend.set_B(type_count,j);
                    lin_bend.set_C(type_count,k);
                    lin_bend.set_linval(type_count,2);
                    lin_bend.set_val(type_count++,0.0);
                  }
                }
            }
      cnt_lin_bend = type_count;

      /* Determine torsions */
      type_count = 0;
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
                      if (cnt) 
                        type_count++;
                      else {
                        tors.set_id(type_count,++id_count);
                        tors.set_A(type_count,a);
                        tors.set_B(type_count,b);
                        tors.set_C(type_count,c);
                        tors.set_D(type_count,d);
                        tors.set_val(type_count++,0.0);
                      }
                    }
                  }
              }
          }

      // add bonus torsions aa-a-b-c-cc if a-b-c is linear, etc.
      int aa, bb, cc, ia, ib, ic, ba, bc, lin_aa, lin_cc, found_aa, found_cc;
      int *taa, *tbb, *tcc, *tdd, skip;
      taa = new int[100]; tbb = new int[100]; tcc = new int[100]; tdd = new int[100];

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
                  fprintf(outfile,"found aa to a,b,c to %d, %d %d %d\n", aa,a,ba,c);


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
                      fprintf(outfile,"found cc to a,b,c to %d %d %d, %d\n", a,bc,c,cc);
                      found_cc = 1;
                      if (is_unique(aa,a,c,cc,taa,tbb,tcc,tdd)) {
                        if (cnt) 
                          type_count++;
                        else {
                          tors.set_id(type_count,++id_count);
                          tors.set_A(type_count,aa);
                          tors.set_B(type_count,a);
                          tors.set_C(type_count,c);
                          tors.set_D(type_count,cc);
                          tors.set_val(type_count++,0.0);
                        }
                      } // include if unique
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
      if ( natom==4 && type_count==0) { 
        int nbonds=0;
        for (a=0;a<natom;++a)
          for (b=0;b<a;++b) {
            if (bonds[a][b])
              nbonds++;
          }
        if (nbonds == 3) {
          // make sure whole molecule is not linear
          if (cnt_bend > 1) {
            // we'll just try 3 torsions for now - this works with NH3 with different
            // atom numberings -- these torsions may have to be changed to
            // use the central atom in a specific way someday
            if (cnt) 
              type_count += 3;
            else {
              tors.set_id(type_count,++id_count);
              a=0;b=1;c=2,d=3;
              swap_tors(&a,&b,&c,&d);
              tors.set_A(type_count,a);
              tors.set_B(type_count,b);
              tors.set_C(type_count,c);
              tors.set_D(type_count,d);
              tors.set_val(type_count++,0.0);
              tors.set_id(type_count,++id_count);

              a=0;b=2;c=3,d=1;
              swap_tors(&a,&b,&c,&d);
              tors.set_A(type_count,a);
              tors.set_B(type_count,b);
              tors.set_C(type_count,c);
              tors.set_D(type_count,d);
              tors.set_val(type_count++,0.0);
              tors.set_id(type_count,++id_count);

              a=0;b=3;c=1,d=2;
              swap_tors(&a,&b,&c,&d);
              tors.set_A(type_count,a);
              tors.set_B(type_count,b);
              tors.set_C(type_count,c);
              tors.set_D(type_count,d);
              tors.set_val(type_count++,0.0);
            }
          } 
        }
      }

      cnt_tors = type_count;
      cnt_out = 0;

      free_block(coord_2d);
      free_block(atom_dist);

        // write out internals newly generated internals
      if (!cnt) {
        ffile(&fp_intco, "intco.dat",0);
        print(fp_intco,0); 
        fclose(fp_intco);
      }
    } // end generate internals loop
    } // end count, creation loop

    num = stre.get_num() + bend.get_num() + tors.get_num() +
      out.get_num() + lin_bend.get_num();
    if (num == 0) {
      fprintf(outfile,"Error: No simple internals were generated and read.  You may ");
      fprintf(outfile,"have to increase the value of the scale_connectivity keyword.");
      exit(2);
    }
    return;
  }



/*-----------------------------------------------------------------------------

  INTERNALS::PRINT

  Print out internal coordinates to a file in intco.dat format
  (print_flag == 0) or with values (print_flag == 1) 
  ----------------------------------------------------------------------------*/

/* Print out internal coordinates to a file in
   intco.dat format (print_flag == 0) or
   with values (print_flag == 1) */
void internals :: print(FILE *fp_out, int print_flag) {
  if (print_flag == 0) fprintf(fp_out, "intco: (\n");
      stre.print(fp_out, print_flag);
      bend.print(fp_out, print_flag);
      tors.print(fp_out, print_flag);
      out.print(fp_out, print_flag);
      lin_bend.print(fp_out, print_flag);
      if (print_flag == 0) fprintf(fp_out, ")\n");
  return;
}



/*----------------------------------------------------------------------------

  INTERNALS::COMPUTE
  ----------------------------------------------------------------------------*/

void internals :: compute_internals(int natom, double *geom) {
  stre.compute(natom, geom);
  bend.compute(natom, geom);
  tors.compute(natom, geom);
  out.compute(natom, geom);
  lin_bend.compute(natom, geom);
  return;
}



/*-----------------------------------------------------------------------------

  INTERNALS::COMPUTE_S
  -----------------------------------------------------------------------------*/

void internals :: compute_s(int natom, double *geom) {
  stre.compute_s(natom, geom);
  bend.compute_s(natom, geom);
  tors.compute_s(natom, geom);
  out.compute_s(natom, geom);
  lin_bend.compute_s(natom, geom);
  return;
}



/*----------------------------------------------------------------------------

  INTERNALS::PRINT_S
  ----------------------------------------------------------------------------*/

void internals :: print_s() {
  stre.print_s();
  bend.print_s();
  tors.print_s();
  out.print_s();
  lin_bend.print_s();
  return;
}


int is_unique(int a, int b, int c, int d, int *aa, int *bb, int *cc, int *dd) {
  int i,unique=1;
  static int n=0;

  for (i=0; i<n; ++i) {
    if ( (a==aa[i]) && (b==bb[i]) && (c==cc[i]) && (d==dd[i]) )
      unique = 0;
  }
  if (unique) {
    aa[n] = a; bb[n] = b; cc[n] = c; dd[n] = d;
    ++n;
  }
  fprintf(outfile,"unique returning %d\n",unique);
  return unique;
}
