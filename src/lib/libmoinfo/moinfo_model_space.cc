#include <iostream>
#include <cmath>
#include <cstdlib>

#include <psifiles.h>
#include <liboptions/liboptions.h>
#include <libutil/libutil.h>

#include "moinfo.h"

extern FILE *outfile;

using namespace std;

namespace psi {

void MOInfo::print_model_space()
{
  fprintf(outfile,"\n");
  fprintf(outfile,"\n  Model space");
  fprintf(outfile,"\n  ------------------------------------------------------------------------------");
  for(int i=0;i<references.size();i++){
    fprintf(outfile,"\n  %2d) ",i);
    references[i].print_occ();
  }
  fprintf(outfile,"\n  ==============================================================================\n");
}

void MOInfo::build_model_space()
{
  /********************************************************
    Generate all the Slater Determinants belonging to the
    model space using the following restrictions:
    -docc are doubly    occupied
    -actv are partially occupied
    -the generalized occupied orbital indexing (docc + actv)
     is assumed (see moinfo.cpp)
  ********************************************************/
  int index;
  MOInfo::SlaterDeterminant docc_det;

  /***********************************************
    Generate all combinations of active orbitals
  ***********************************************/
  if(nactv_docc == 0){
    /********************************************************
      Set up the doubly occupied part of all the
      determinants in the model space
    ********************************************************/
    index = 0;
    for(int h = 0; h < nirreps; ++h){
      for(int i = 0; i < docc[h]; ++i){
        docc_det.set(index);
        docc_det.set(index + nall);
        index++;
      }
      index += actv[h];
      index += extr[h];
    }

    /********************************************************
      Set up the a vectors containing the active orbitals and
      their symmetry
    ********************************************************/
    std::vector<int> alpha_active,alpha_active_sym,beta_active,beta_active_sym;
    index = 0;
    for(int h = 0; h < nirreps; ++h){
      index += docc[h];
      for(int i = 0; i < actv[h]; ++i){
        alpha_active.push_back(index);
        alpha_active_sym.push_back(h);
        beta_active.push_back(index + nall);
        beta_active_sym.push_back(h);
        index++;
      }
      index += extr[h];
    }

    std::vector<std::vector<int> > alpha_combinations,beta_combinations;
    generate_combinations(nactv,nactive_ael,alpha_combinations);
    generate_combinations(nactv,nactive_bel,beta_combinations);
    if(alpha_combinations.size()==0)
      alpha_combinations.push_back(vector<int>(0));
    if(beta_combinations.size()==0)
      beta_combinations.push_back(vector<int>(0));
    for(int a=0;a<alpha_combinations.size();a++){
      for(int b=0;b<beta_combinations.size();b++){
        int sym = 0; // Symmetry of the determinant
        // Create a copy of the docc_det
        SlaterDeterminant det(docc_det);
        // Fill the alpha active orbitals
        for(int i=0;i<nactive_ael;i++){
          det.set(alpha_active[alpha_combinations[a][i]]);
          sym = sym ^ alpha_active_sym[alpha_combinations[a][i]];
        }
        // Fill the beta active orbitals
        for(int i=0;i<nactive_bel;i++){
          det.set(beta_active[beta_combinations[b][i]]);
          sym = sym ^  beta_active_sym[beta_combinations[b][i]];
        }
        // Test the wfn symmetry
        if(sym==wfn_sym){
          if(det.is_closed_shell()){
            // Closed-shell determinant
            closed_shell_refs.push_back(references.size());
            unique_refs.push_back(references.size());
            all_refs.push_back(references.size());
          }else{
            // Open-shell determinant
            bool add_it = true;
            int  spin_mirror = references.size();
            if(options_get_bool("USE_SPIN_SYMMETRY")){
              // Check if this is a spin-flipped determinant
              for(int ref=0;ref<references.size();ref++){
                if(references[ref].is_spin_flipped(det)){
                  add_it      = false;
                  spin_mirror = ref;
                }
              }
              if(add_it){
                unique_open_shell_refs.push_back(references.size());
                unique_refs.push_back(references.size());
                all_refs.push_back(references.size());
              }else{
                all_refs.push_back(spin_mirror);
              }
            }
          }
          references.push_back(det);
        }
      }
    }
  }else{
    /********************************************************
      Set up the doubly occupied part of all the
      determinants in the model space
    ********************************************************/
    index = 0;
    for(int h=0;h<nirreps;h++){
      for(int i=0;i<docc[h] + actv_docc[h];i++){
        docc_det.set(index);
        docc_det.set(index + nall);
        index++;
      }
      index += actv[h] - actv_docc[h];
      index += extr[h];
    }
    closed_shell_refs.push_back(references.size());
    unique_refs.push_back(references.size());
    all_refs.push_back(references.size());
    references.push_back(docc_det);
  }

  if(references.size() == 0){
    fprintf(outfile,"\n\n  MOInfo found no reference in the model space");
    fprintf(outfile,"\n  Please check the following:");
    fprintf(outfile,"\n  1) Definition of FOCC, DOCC, ACTV, and FVIR");
    fprintf(outfile,"\n  2) Symmetry of the wavefunction");
    fprintf(outfile,"\n  3) Charge and multiplicity");
    fprintf(outfile,"\n\n  PSIMRCC will end the computation.\n");
    fflush(outfile);
    exit(PSI_RETURN_FAILURE);
  }
}

/*!
    \fn MOInfo::make_internal_excitations()
 */
void MOInfo::make_internal_excitations()
{
  /****************************************
    Build the mappings between references
    |m> = (+/-) ... b+ j a+ i |n>
  ****************************************/
  for(int m=0;m<references.size();m++){
    vector<vector<pair<int,int> > > alpha_internals_ref_m;
    vector<vector<pair<int,int> > >  beta_internals_ref_m;
    vector<double>                   sign_internals_ref_m;
//     DEBUGGING(1,
//       fprintf(outfile,"\n\n\tReference ");
//       references[m].print_occ();
//       fprintf(outfile," gives:");
//     );
    for(int n=0;n<references.size();n++){
      double sign=1.0;
      std::vector<pair<int,int> > alpha_operators;
      std::vector<pair<int,int> > beta_operators;
      references[m].get_internal_excitations(references[n],sign,alpha_operators,beta_operators);
      alpha_internals_ref_m.push_back(alpha_operators);
      beta_internals_ref_m.push_back(beta_operators);
      sign_internals_ref_m.push_back(sign);
//       DEBUGGING(1,
//         fprintf(outfile,"\n\t  ");
//         references[n].print_occ();
//         fprintf(outfile," = %s{",sign > 0.0 ? "+" : (sign == 0.0 ? "0" : "-"));
//         for(int i = 0; i<beta_operators.size();i++)
//           fprintf(outfile," %db+ %db-",beta_operators[i].second,beta_operators[i].first);
//         for(int i = 0; i<alpha_operators.size();i++)
//           fprintf(outfile," %da+ %da-",alpha_operators[i].second,alpha_operators[i].first);
//         fprintf(outfile," }");
//         references[m].print_occ();
//       );
    }
    alpha_internal_excitations.push_back(alpha_internals_ref_m);
    beta_internal_excitations.push_back(beta_internals_ref_m);
    sign_internal_excitations.push_back(sign_internals_ref_m);
  }
}

vector<int> MOInfo::get_aocc(string str,int i)
{
  int i_ref = get_ref_number(str,i);
  return(references[i_ref].get_aocc());
}

vector<int> MOInfo::get_bocc(std::string str,int i)
{
  int i_ref = get_ref_number(str,i);
  return(references[i_ref].get_bocc());
}

vector<int> MOInfo::get_avir(std::string str,int i)
{
  int i_ref = get_ref_number(str,i);
  return(references[i_ref].get_avir());
}

vector<int> MOInfo::get_bvir(std::string str,int i)
{
  int i_ref = get_ref_number(str,i);
  return(references[i_ref].get_bvir());
}

vector<int> MOInfo::get_aocc(int i)
{
  return(references[i].get_aocc());
}

vector<int> MOInfo::get_bocc(int i)
{
  return(references[i].get_bocc());
}

vector<int> MOInfo::get_auoc(int i)
{
  return(references[i].get_avir());
}

vector<int> MOInfo::get_buoc(int i)
{
  return(references[i].get_bvir());
}

vector<int> MOInfo::get_determinant(int i)
{
  vector<int> occupation(nall * 2,0);
  for(int p = 0; p < 2 * nall; ++p)
    if(references[i].test(p))
      occupation[p] = 1;
  return occupation;
}

vector<pair<int,int> > MOInfo::get_alpha_internal_excitation(int i,int j)
{
  return(alpha_internal_excitations[i][j]);
}

vector<pair<int,int> > MOInfo::get_beta_internal_excitation(int i,int j)
{
  return(beta_internal_excitations[i][j]);
}

double  MOInfo::get_sign_internal_excitation(int i,int j)
{
  return(sign_internal_excitations[i][j]);
}

/*!
    \fn MOInfo::get_ref_number(string str, int n)
 */
int MOInfo::get_ref_number(string str, int n)
{
  if(str=="a")
    return(all_refs[n]);
  if(str=="u")
    return(unique_refs[n]);
  if(str=="c")
    return(closed_shell_refs[n]);
  if(str=="o")
    return(unique_open_shell_refs[n]);
  print_error(outfile,"MOInfo::get_ref_number(string str, int n) undefined space", __FILE__,__LINE__);
  return(NULL);
}

/*!
    \fn MOInfo::get_ref_size(string str)
 */
int MOInfo::get_ref_size(string str)
{
  if(str=="a")
    return(all_refs.size());
  if(str=="u")
    return(unique_refs.size());
  if(str=="c")
    return(closed_shell_refs.size());
  if(str=="o")
    return(unique_open_shell_refs.size());
  print_error(outfile,"MOInfo::get_ref_size(string str) undefined space", __FILE__,__LINE__);
  return(NULL);
}

vector<string> MOInfo::get_matrix_names(std::string str)
{
  vector<string> names;
  if(str.find("{a}")!=string::npos){
    for(int n=0;n<all_refs.size();n++)
      names.push_back(find_and_replace(str,"{a}","{" + to_string(all_refs[n]) +"}"));
  }else if(str.find("{u}")!=string::npos){
    for(int n=0;n<unique_refs.size();n++)
      names.push_back(find_and_replace(str,"{u}","{" + to_string(unique_refs[n]) +"}"));
  }else if(str.find("{c}")!=string::npos){
    for(int n=0;n<closed_shell_refs.size();n++)
      names.push_back(find_and_replace(str,"{c}","{" + to_string(closed_shell_refs[n]) +"}"));
  }else if(str.find("{o}")!=string::npos){
    for(int n=0;n<unique_open_shell_refs.size();n++)
      names.push_back(find_and_replace(str,"{o}","{" + to_string(unique_open_shell_refs[n]) +"}"));
  }else
    names.push_back(str);
  return(names);
}

}
