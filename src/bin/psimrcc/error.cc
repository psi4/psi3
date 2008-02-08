/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/



#ifndef IOSTREAM_H
#define IOSTREAM_H
#include <iostream>
#endif
#include "error.h"

extern FILE *infile, *outfile;

namespace psi{ namespace psimrcc{

using namespace std;

void print_error(string message, char* file, int line)
{
  fprintf(outfile,"\n\n%s in file %s, line %d\n",message.c_str(),file,line);
  fflush(outfile);
  exit(1);
}

void print_error(char* message, char* file, int line)
{
  print_error(message,file,line,1);
}

void print_error(char* message, char* file, int line, int error)
{
  fprintf(outfile,"\n\n%s in file %s, line %d\n",message,file,line);
  fflush(outfile);
  exit(error);
}

void print_developing(char* message, char* file, int line)
{
  print_developing(message,file,line,1);
}

void print_developing(char* message, char* file, int line, int error)
{
  fprintf(outfile,"\n\n%s: feature not yet implemented in file %s, line %d\n",message,file,line);
  fflush(outfile);
  exit(error);
}

}} /* End Namespaces */