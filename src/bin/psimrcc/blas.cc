#include <liboptions/liboptions.h>
#include "blas.h"
#include "memory_manager.h"
#include "debugging.h"
#include <libmoinfo/libmoinfo.h>
#include <libutil/libutil.h>
#include <algorithm>

#include <libciomr/libciomr.h>
#include <libqt/qt.h>

extern FILE *infile, *outfile;

namespace psi{ namespace psimrcc{

using namespace std;

CCBLAS::CCBLAS():
full_in_core(false),work_size(0)
{
  init();
}


CCBLAS::~CCBLAS()
{
  cleanup();
}

void CCBLAS::init()
{
  add_indices();
  allocate_work();
  allocate_buffer();
}

void CCBLAS::cleanup()
{
  free_sortmap();
  free_buffer();
  free_work();
  free_matrices();
  free_indices();
}

void CCBLAS::allocate_work()
{
  // Make sure work is empty
  if(!work.empty())
    for(int n=0;n<work.size();n++)
      if(work[n]!=NULL)
        delete[] work[n];

  for(int n=0;n<options_get_int("NUM_THREADS");n++)
    work.push_back(NULL);
  // Compute the temporary work space size
  CCIndex* oo_pair = get_index("[oo]");
  CCIndex* vv_pair = get_index("[vv]");
  work_size = 0;
  for(int h=0;h<moinfo->get_nirreps();h++){
    if(oo_pair->get_pairpi(h)<=vv_pair->get_pairpi(h))
      work_size += oo_pair->get_pairpi(h)*vv_pair->get_pairpi(h);
    else
      work_size += oo_pair->get_pairpi(h)*oo_pair->get_pairpi(h);
  }
  // Allocate the temporary work space
  for(int n=0;n<options_get_int("NUM_THREADS");n++){
    allocate1(double,work[n],work_size);
//     work[n]=new double[work_size]; MEM_
    zero_arr(work[n],work_size);
  }
  mem->add_allocated_memory(to_MB(work_size));
}

void CCBLAS::allocate_buffer()
{
  // Make sure buffer is empty
  if(!buffer.empty())
    for(int n=0;n<buffer.size();n++)
      if(buffer[n]!=NULL)
        delete[] buffer[n];

  for(int n=0;n<options_get_int("NUM_THREADS");n++)
    buffer.push_back(NULL);
  // Compute the temporary buffer space size, 101% of the actual strip size
  buffer_size = 1.01*mem->get_integral_strip_size() / to_MB(1);
  // Allocate the temporary buffer space
  for(int n=0;n<options_get_int("NUM_THREADS");n++){
    buffer[n]=new double[buffer_size];
    zero_arr(buffer[n],buffer_size);
  }
  mem->add_allocated_memory(to_MB(buffer_size));
}

void CCBLAS::free_sortmap()
{
  for(SortMap::iterator iter=sortmap.begin();iter!=sortmap.end();++iter){
    for(int irrep=0;irrep<moinfo->get_nirreps();irrep++)
      delete[] iter->second[irrep];
    delete[] iter->second;
  }
}

void CCBLAS::free_work()
{
  // Delete the temporary work space
  for(int n=0;n<work.size();n++){
    if(work[n]!=NULL){
      release1(work[n]);
//        delete[] work[n]; MEM_
    }
  }
}

void CCBLAS::free_buffer()
{
  // Delete the temporary buffer space
  for(int n=0;n<buffer.size();n++){
    if(buffer[n]!=NULL){
      delete[] buffer[n];
    }
  }
}

void CCBLAS::free_matrices()
{
  for(MatrixMap::iterator iter=matrices.begin();iter!=matrices.end();++iter){
    delete iter->second;
  }
}

void CCBLAS::free_indices()
{
  for(IndexMap::iterator iter=indices.begin();iter!=indices.end();++iter){ 
    delete iter->second;
  }
}

void CCBLAS::add_indices()
{
  add_index("[]");
  add_index("[o]");
  add_index("[v]");
  add_index("[a]");
  add_index("[o>o]");
  add_index("[v>v]");
  add_index("[v>=v]");
  add_index("[oo]");
  add_index("[ov]");
  add_index("[vo]");
  add_index("[vv]");
  add_index("[aa]");
  add_index("[aaa]");
  add_index("[ooo]");
  add_index("[oov]");
  add_index("[voo]");
  add_index("[ovv]");
  add_index("[vvo]");
  add_index("[ovo]");

  // MP2-CCSD
  if(options_get_str("CORR_WFN")=="MP2-CCSD"){
    add_index("[oav]");
    add_index("[ova]");
    add_index("[avo]");
    add_index("[aao]");
    add_index("[aoa]");
    add_index("[oaa]");
    add_index("[vaa]");
    add_index("[aav]");
    add_index("[ava]");
  }
  if(options_get_str("CORR_WFN")!="PT2"){
    add_index("[vvv]");
  }

  // Mk-MRPT2
  add_index("[ao]");
  add_index("[av]");

  // Not useful
  add_index("[oa]");
  add_index("[va]");
}

// void CCBLAS::allocate_matrices_in_core()
// {
//   for(MatrixMap::iterator iter=matrices.begin();iter!=matrices.end();++iter){
//     CCMatrix* Matrix = iter->second;
//     fprintf(outfile,"\n%s(analyzing)",Matrix->get_label().c_str());
//     fflush(outfile);
//     if(Matrix->get_out_of_core()){
//       Matrix->load();
//       fprintf(outfile,"\n%s <- reading from disk",Matrix->get_label().c_str());
//       fflush(outfile);
//     }else if(!Matrix->is_allocated())
//       Matrix->allocate_memory();
//   }
// }

void CCBLAS::print(const char* cstr)
{
  string str(cstr);
  vector<string> names = moinfo->get_matrix_names(str);
  for(int n=0;n<names.size();n++)
    print_ref(names[n]);
}

void CCBLAS::print_ref(string& str)
{
  get_Matrix(str)->print();
}

void CCBLAS::print_memory()
{
  double total_memory_required =0.0;
  fprintf(outfile,"\n\n\t-----------------------------------------------------------------------------");
  fprintf(outfile,"\n\tMatrix ID    Memory(MB)   Cumulative Memory(MB)  Accessed    Label");
  fprintf(outfile,"\n\t------------------------------------------------------------------------------");

  for(MatrixMap::iterator iter=matrices.begin();iter!=matrices.end();++iter){
    total_memory_required += iter->second->get_memory();
    fprintf(outfile,"\n\t%4d",distance(matrices.begin(),iter));
    fprintf(outfile,"     %10.2f",iter->second->get_memory());
    fprintf(outfile,"        %10.2f",total_memory_required);
    fprintf(outfile,"             %4d",iter->second->get_naccess());
    fprintf(outfile,"         %s",iter->second->get_label().c_str());
  }
  fprintf(outfile,"\n\t------------------------------------------------------------------------------");
  fprintf(outfile,"\n\n\tTotal memory required for matrices = %10.2f (Mb)\n",total_memory_required);

  total_memory_required =0.0;

  fprintf(outfile,"\n\n\t-------------------------------------------------------------");
  fprintf(outfile,"\n\tIndex ID    Memory(MB)   Cumulative Memory(MB)     Label");
  fprintf(outfile,"\n\t--------------------------------------------------------------");

  for(IndexMap::iterator iter=indices.begin();iter!=indices.end();++iter){
    total_memory_required += iter->second->get_memory();
    fprintf(outfile,"\n\t%4d",distance(indices.begin(),iter));
    fprintf(outfile,"     %10.2f",iter->second->get_memory());
    fprintf(outfile,"        %10.2f",total_memory_required);
    fprintf(outfile,"         %s",iter->second->get_label().c_str());
  }
  fprintf(outfile,"\n\t--------------------------------------------------------------");

  fprintf(outfile,"\n\n\tTotal memory required for indexing = %10.2f (Mb)\n",total_memory_required);
}

/**
 * This routine computes which quantities have to be initially stored in memory and which on disk
 */
int CCBLAS::compute_storage_strategy()
{
  fprintf(outfile,"\n#CC ----------------------------------");
  fprintf(outfile,"\n#CC    Computing Storage Strategy");
  fprintf(outfile,"\n#CC ----------------------------------");

  double available_memory     = mem->get_free_memory();
  double fraction_for_in_core = 0.97; // Fraction of the total available memory that may be used
  double storage_memory       = available_memory * fraction_for_in_core;
  double fully_in_core_memory = 0.0;
  double integrals_memory     = 0.0;
  double fock_memory          = 0.0;
  double others_memory        = 0.0;

  fprintf(outfile,"\n#CC Input memory                                = %10.2f Mb",mem->get_total_memory());
  fprintf(outfile,"\n#CC Free memory                                 = %10.2f Mb",available_memory);
  fprintf(outfile,"\n#CC Free memory available for matrices          = %10.2f Mb (%3.0f%s of free_memory)",storage_memory,fraction_for_in_core*100.0,"%");

  // Gather the memory requirements for all the CCMAtrix object
  // and divide the integrals from all the other matrices.
  // At the same time compute the memory requirements for
  // a fully in-core algorithm.
  vector<pair<double,pair<CCMatrix*,int> > > integrals;
  vector<pair<double,pair<CCMatrix*,int> > > fock;
  vector<pair<double,pair<CCMatrix*,int> > > others;
  for(MatrixMap::iterator it=matrices.begin();it!=matrices.end();++it){
    for(int h=0;h<moinfo->get_nirreps();++h){
      double block_memory = it->second->get_memorypi(h);
      if(it->second->is_integral()){
        integrals.push_back(make_pair(block_memory,make_pair(it->second,h)));
        integrals_memory += block_memory;
      }else if(it->second->is_fock()){
        fock.push_back(make_pair(block_memory,make_pair(it->second,h)));
        fock_memory += block_memory;
      }else{
        others.push_back(make_pair(block_memory,make_pair(it->second,h)));
        others_memory += block_memory;
      }
      fully_in_core_memory += block_memory;
    }
  }
  fprintf(outfile,"\n#CC Memory required by fock matrices            = %10.2f Mb",fock_memory);
  fprintf(outfile,"\n#CC Memory required by integrals                = %10.2f Mb",integrals_memory);
  fprintf(outfile,"\n#CC Memory required by other matrices           = %10.2f Mb",others_memory);
  fprintf(outfile,"\n#CC Memory required for full in-core algorithm  = %10.2f Mb",fully_in_core_memory);

  // Check if you may use a fully in core algorithm
  full_in_core = false;
  int strategy = 0;
  if(fully_in_core_memory < storage_memory ){
    full_in_core = true;
    fprintf(outfile,"\n#CC PSIMRCC will perform a full in-core computation");
    strategy = 0;
  }else{
    if(others_memory < storage_memory ){
      fprintf(outfile,"\n#CC PSIMRCC will store some integrals out-of-core");
      strategy = 1;
    }else{
      fprintf(outfile,"\n#CC PSIMRCC will store all integrals and some other matrices out-of-core");
      strategy = 2;
      fprintf(outfile,"\n#CC CCBLAS::compute_storage_strategy(): Strategy #2 is not implemented yet");
      fflush(outfile);
      exit(1);
    }
  }
  sort(integrals.begin(),integrals.end());
  sort(others.begin(),others.end());
  for(int i=0;i<fock.size();i++){
    // Store all the fock matrices in core and allocate them
    storage_memory -= fock[i].first;
    load_irrep(fock[i].second.first,fock[i].second.second);
  }
  // Let the CCBlas class worry about allocating matrices
  int number_of_others_on_disk = 0;
  for(int i=0;i<others.size();i++){
    // Check if this matrix can be stored in core
    if(others[i].first < storage_memory){
      storage_memory -= others[i].first;
      load_irrep(others[i].second.first,others[i].second.second);
    }else{
      number_of_others_on_disk++;
    }
  }
  int number_of_integrals_on_disk = 0;
  for(int i=0;i<integrals.size();i++){
    // Check if this matrix can be stored in core
    if(integrals[i].first < storage_memory){
      storage_memory -= integrals[i].first;
      load_irrep(integrals[i].second.first,integrals[i].second.second);
    }else{
      number_of_integrals_on_disk++;
    }
  }

  DEBUGGING(1,
    fprintf(outfile,"\n#CC -------------------- Fock matrices -------------------------");
    for(int i=0;i<fock.size();i++){
      if(fock[i].first > 1.0e-5){
        fprintf(outfile,"\n#CC %-32s irrep %d   %6.2f Mb --> ",fock[i].second.first->get_label().c_str(),
                                                      fock[i].second.second,
                                                      fock[i].first);
        fprintf(outfile,"%s",fock[i].second.first->is_block_allocated(fock[i].second.second) ? "in-core" : "out-of-core");
      }
    }
    fprintf(outfile,"\n#CC -------------------- Other matrices ------------------------");    
    for(int i=0;i<others.size();i++){
      if(others[i].first > 1.0e-5){
        fprintf(outfile,"\n#CC %-32s irrep %d   %6.2f Mb --> ",others[i].second.first->get_label().c_str(),
                                                      others[i].second.second,
                                                      others[i].first);
        fprintf(outfile,"%s",others[i].second.first->is_block_allocated(others[i].second.second) ? "in-core" : "out-of-core");
      }
    }
    fprintf(outfile,"\n#CC -------------------- Integrals -----------------------------");
    for(int i=0;i<integrals.size();i++){
      if(integrals[i].first > 1.0e-5){
        fprintf(outfile,"\n#CC %-32s irrep %d   %6.2f Mb --> ",integrals[i].second.first->get_label().c_str(),
                                                      integrals[i].second.second,
                                                      integrals[i].first);
        fprintf(outfile,"%s",integrals[i].second.first->is_block_allocated(integrals[i].second.second) ? "in-core" : "out-of-core");
      }
    }
    fprintf(outfile,"\n\n");
  );

  if(!full_in_core){
    fprintf(outfile,"\n#CC Out-of-core algorithm will store %d other matrices on disk",number_of_others_on_disk);
    fprintf(outfile,"\n#CC Out-of-core algorithm will store %d integrals on disk",number_of_integrals_on_disk);
  }
  return(strategy);
}

/**
 * This routine computes which quantities have to be initially stored in memory and which on disk
 */
void CCBLAS::show_storage()
{
  DEBUGGING(1,
    fprintf(outfile,"\n#CC ----------------------------------");
    fprintf(outfile,"\n#CC            Show Storage ");
    fprintf(outfile,"\n#CC ----------------------------------");

    for(MatrixMap::iterator it=matrices.begin();it!=matrices.end();++it){
      for(int h=0;h<moinfo->get_nirreps();++h){
        double block_memory = it->second->get_memorypi(h);
        fprintf(outfile,"\n#CC %-32s irrep %d   %6.2f Mb",it->second->get_label().c_str(),h,block_memory);
        fprintf(outfile," is %s",it->second->is_block_allocated(h) ? "allocated" : "not allocated");
        fprintf(outfile,"%s",it->second->is_out_of_core(h) ? "(out-of-core)" : "");
  
      }
    }
  )
}

}} /* End Namespaces */
