#include <cstring>
#include <iostream>

#include <libqt/qt.h>

#include "vector_base.h"
#include "memory_manager.h"

extern FILE* outfile;

VectorBase::VectorBase(int elements) : elements_(elements),vector_(NULL)
{
  allocate1(double,vector_,elements_);
}

VectorBase::~VectorBase()
{
  release1(vector_);
}

void VectorBase::print()
{
  fprintf(outfile,"\n  ");
  for(int i=0 ; i < elements_; i++){
    fprintf(outfile,"%10.6f",vector_[i]);
  }
}
