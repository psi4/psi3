#ifndef _psi_src_bin_psimrcc_error_h
#define _psi_src_bin_psimrcc_error_h
/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

#include <string>

namespace psi{ namespace psimrcc{

void print_error(std::string message, char* file, int line);
void print_error(char* message, char* file, int line);
void print_error(char* message, char* file, int line,int error);
void print_developing(char* message, char* file, int line);
void print_developing(char* message, char* file, int line,int error);

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_error_h