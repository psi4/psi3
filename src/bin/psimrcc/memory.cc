/***************************************************************************
 *   Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *   frank@ccc.uga.edu
 *   SR/MRCC Code
 ***************************************************************************/
// #include "ccmemory.h"
// 
// CCMemoryManager::CCMemoryManager()
// {
// }
// 
// CCMemoryManager::~CCMemoryManager()
// {
//   free_memory();
// }
// 
// CCMemoryBlock<double*>& CCMemoryManager::init_block(size_t size1)
// {
//   // Allocate the memory
//   std::vector<size_t> size;
//   size.push_back(size1);
//   size_t size_product = size1;
//   double* block = new double[size_product];
//   one_index.push_back(CCMemoryBlock<double*>(block,size));
//   return(one_index[one_index.size()]);
// }
// 
// void CCMemoryManager::free_memory()
// {
//   OneVecIt one_end_it = one_index.end();
//   for(OneVecIt one_it = one_index.begin();one_it!=one_end_it;++one_it){
//     
//   }
// }