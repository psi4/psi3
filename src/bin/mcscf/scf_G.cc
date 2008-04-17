#include <iostream>

#include <psifiles.h>
#include <libpsio/psio.h>

#include "scf.h"
#include "memory_manager.h"

extern FILE* outfile;

namespace psi{ namespace mcscf{

void SCF::construct_G(SBlockMatrix& density,
                      SBlockMatrix& G,
                      double* integrals,
                      int batch)
{
  construct_G(density,G,integrals,batch,1.0);
}


void SCF::construct_G(SBlockMatrix& density,
                      SBlockMatrix& G,
                      double* integrals,
                      int batch,
                      double factor)
{
  double* D_vector;
  double* G_vector;

  allocate1(double,D_vector,npairs);
  allocate1(double,G_vector,npairs);

  // Convert D to a vector and double the off diagonal elements
  for(int h = 0; h < nirreps; ++h){
    for(int p = 0; p < sopi[h]; ++p){
      int p_abs = p + block_offset[h];
      for(int q = 0; q <= p; ++q){
        int q_abs = q + block_offset[h];
        D_vector[ pair[p_abs][q_abs] ] = 2.0 * density->get(h,p,q);
        G_vector[ pair[p_abs][q_abs] ] = 0.0;
      }
      D_vector[ pair[p_abs][p_abs] ]  *= 0.5;
    }
  }

  // general algorithm
  double G_pq,D_pq;
  double* D_rs;
  double* G_rs;
  int     pq,rs;
  double* PK_pqrs = integrals;

  for(pq = batch_pq_min[batch]; pq < batch_pq_max[batch]; ++pq){
    G_pq = 0.0;
    D_pq = D_vector[pq];
    D_rs = &D_vector[0];
    G_rs = &G_vector[0];
    for(rs = 0; rs <= pq; ++rs){
      G_pq += *PK_pqrs * (*D_rs);
      *G_rs += *PK_pqrs * D_pq;
      ++D_rs;
      ++G_rs;
      ++PK_pqrs;
    }
    G_vector[pq] += G_pq;
  }

  // Convert G to a matrix
  for(int h = 0; h < nirreps; ++h){
    for(int p = 0; p < sopi[h]; ++p){
      int p_abs = p + block_offset[h];
      for(int q = 0; q < sopi[h]; ++q){
        int q_abs = q + block_offset[h];
        G->set(h,p,q,2.0 * factor * G_vector[ pair[p_abs][q_abs]]);
      }
    }
  }

  release1(G_vector);
  release1(D_vector);
}



void SCF::read_Raffanetti(const char* integral_type, double* integrals, int batch)
{
  // Read the PK matrix from disk
  char data_label[80];
  sprintf(data_label,"%s_%d",integral_type,batch);
  size_t buffer_size = batch_size[batch] * sizeof(double);
  psio_read_entry(PSIF_MCSCF,data_label,(char*)integrals,buffer_size);
}

void SCF::write_Raffanetti(const char* integral_type, double* integrals, int batch)
{
  // Write the PK matrix to disk
  char data_label[80];
  sprintf(data_label,"%s_%d",integral_type,batch);
  size_t buffer_size = batch_size[batch] * sizeof(double);
  psio_write_entry(PSIF_MCSCF,data_label,(char*)integrals,buffer_size);
}


}} /* End Namespaces */


/*
void SCF::construct_G(SBlockMatrix& density, SBlockMatrix& G, const char* integral_type, double factor)
{
  double* D_vector;
  double* G_vector;

  allocate1(double,D_vector,npairs);
  allocate1(double,G_vector,npairs);

  // Convert D to a vector and double the off diagonal elements
  for(int h = 0; h < nirreps; ++h){
    for(int p = 0; p < sopi[h]; ++p){
      int p_abs = p + block_offset[h];
      for(int q = 0; q <= p; ++q){
        int q_abs = q + block_offset[h];
        D_vector[ pair[p_abs][q_abs] ] = 2.0 * density->get(h,p,q);
        G_vector[ pair[p_abs][q_abs] ] = 0.0;
      }
      D_vector[ pair[p_abs][p_abs] ] *= 0.5;
    }
  }

  // general algorithm
  double G_pq,D_pq;
  double* D_rs;
  double* G_rs;
  int     pq,rs;

  for(int batch = 0; batch < nbatch; ++batch){
    size_t min_index   = batch_index_min[batch];
    size_t max_index   = batch_index_max[batch];
    size_t buffer_size = max_index - min_index;
    // Read the PK matrix from disk
    char data_label[80];
    sprintf(data_label,"%s_%d",integral_type,batch);
    psio_read_entry(PSIF_MCSCF,data_label,(char*)&(PK[0]),buffer_size * sizeof(double));
    double* PK_pqrs = &(PK[0]);

    for(pq = batch_pq_min[batch]; pq < batch_pq_max[batch]; ++pq){
      G_pq = 0.0;
      D_pq = D_vector[pq];
      D_rs = &D_vector[0];
      G_rs = &G_vector[0];
      for(rs = 0; rs <= pq; ++rs){
        G_pq += *PK_pqrs * (*D_rs);
        *G_rs += *PK_pqrs * D_pq;
        ++D_rs;
        ++G_rs;
        ++PK_pqrs;
      }
      G_vector[pq] += G_pq;
    }
  }

  // Convert G to a matrix
  for(int h = 0; h < nirreps; ++h){
    for(int p = 0; p < sopi[h]; ++p){
      int p_abs = p + block_offset[h];
      for(int q = 0; q < sopi[h]; ++q){
        int q_abs = q + block_offset[h];
        G->set(h,p,q,2.0 * factor * G_vector[ pair[p_abs][q_abs]]);
      }
    }
  }

  release1(G_vector);
  release1(D_vector);
}*/

/////////////////////////////////////////////////
// Routines for multithreaded code
/////////////////////////////////////////////////
// typedef struct{
//   int       thread;
//   int       nthreads;
//   int       total_symmetric_pairpi;
//   pthread_t thread_id;
//   double*   D_vector;
//   double*   G_vector;
//   double*   PK;
// } PK_arguments;
// 
// void* construct_G_PK_thread(void* arguments);
// void multithreaded_construct_G_PK();
// void SCF::multithreaded_construct_G_PK()
// {
//   double* D_vector;
//   allocate1(double,D_vector,npairs);
// 
//   // Convert D to a vector and double the off diagonal elements
//   for(int h = 0; h < nirreps; ++h){
//     for(int p = 0; p < sopi[h]; ++p){
//       int p_abs = p + block_offset[h];
//       for(int q = 0; q <= p; ++q){
//         int q_abs = q + block_offset[h];
//         D_vector[ pair[p_abs][q_abs] ] = 2.0 * D->get(h,p,q);
//       }
//       D_vector[ pair[p_abs][p_abs] ] *= 0.5;
//     }
//   }
// 
//   PK_arguments* PK_args = new PK_arguments[nthreads];
// 
//   // Allocate a G_vector for each thread
//   // Setup PK_args
//   // Start the thread
//   double** G_vector_threads;
//   allocate1(double*,G_vector_threads,nthreads);
// 
//   for(int n = 0 ; n< nthreads ; ++n){
//     allocate1(double,G_vector_threads[n],npairs);
//     PK_args[n].thread   = n;
//     PK_args[n].nthreads = nthreads;
//     PK_args[n].total_symmetric_pairpi = pairpi[0];
//     PK_args[n].D_vector = D_vector;
//     PK_args[n].G_vector = G_vector_threads[n];
//     PK_args[n].PK       = PK;
// 
//     pthread_create(&PK_args[n].thread_id,
//                    NULL,
//                    &construct_G_PK_thread,
//                    &PK_args[n]);
//   }
// 
//   // Join the threads
//   for(int n = 0 ; n< nthreads ; ++n){
//     pthread_join(PK_args[n].thread_id,NULL);
//   }
// 
//   // Sum up all the contributions from the G to a matrices
//   G->zero();
//   for(int n = 0 ; n< nthreads ; ++n){
//     for(int h = 0; h < nirreps; ++h){
//       for(int p = 0; p < sopi[h]; ++p){
//         int p_abs = p + block_offset[h];
//         for(int q = 0; q < sopi[h]; ++q){
//           int q_abs = q + block_offset[h];
//           G->add(h,p,q,2.0 * G_vector_threads[n][ pair[p_abs][q_abs] ] );
//         }
//       }
//     }
//     release1(G_vector_threads[n]);
//   }
// 
//   release1(G_vector_threads);
//   release1(D_vector);
//   delete[] PK_args;
// }
// 
// void* construct_G_PK_thread(void* arguments)
// {
//   PK_arguments* PK_args = static_cast<PK_arguments*>(arguments);
//   int     thread   = PK_args->thread;
//   int     nthreads = PK_args->nthreads;
//   int     total_symmetric_pairpi = PK_args->total_symmetric_pairpi;
//   double  G_pq,D_pq;
//   double* PK_block = PK_args->PK;
//   double* D_vector = PK_args->D_vector;
//   double* G_vector = PK_args->G_vector;
//   double* D_rs;
//   double* G_rs;
//   int pq,rs;
//   for(pq = 0; pq < total_symmetric_pairpi; ++pq){
//     if(pq % nthreads == thread){
//       G_pq = 0.0;
//       D_pq = D_vector[pq];
//       D_rs = &D_vector[0];
//       G_rs = &G_vector[0];
//       for(rs = 0; rs <= pq; ++rs){
//         G_pq += *PK_block * (*D_rs);
//         *G_rs += *PK_block * D_pq;
//         ++D_rs;
//         ++G_rs;
//         ++PK_block;
//       }
//       G_vector[pq] += G_pq;
//     }else{
//       PK_block += pq + 1;
//     }
//   }
//   return(NULL);
// }
///////////////////////////////////////////////////
