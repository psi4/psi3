#ifndef _psi_src_bin_psimrcc_ccindex_iterator_h
#define _psi_src_bin_psimrcc_ccindex_iterator_h

/*! \file    index_iterator.h
    \ingroup (PSIMRCC)
    \brief   This class is used to iterate over n-tuples of MOs indices (p,q,r,..)
*/

#include <string>

namespace psi{ namespace psimrcc{

class CCIndex;

class CCIndexIterator{
public:
  ///////////////////////////////////////////////////////////////////////////////
  // Class Constructor and Destructor
  ///////////////////////////////////////////////////////////////////////////////
  CCIndexIterator(std::string str,int select_irrep = -1);
  CCIndexIterator(CCIndex* index,int select_irrep = -1);
  ~CCIndexIterator();
  ///////////////////////////////////////////////////////////////////////////////
  // Class Methods
  ///////////////////////////////////////////////////////////////////////////////
  bool operator++();
  bool next_irrep();
  bool next_element_in_irrep();

  void        reset();
  void        set_irrep(int n)  {select = n;}

  ///////////////////////////////////////////////////////////////////////////////
  // Class Public Members
  ///////////////////////////////////////////////////////////////////////////////
  int                               sym;
  size_t                            rel;
  size_t                            abs;
  short*                            ind_abs;
  int*                              ind_sym;
private:
  ///////////////////////////////////////////////////////////////////////////////
  // Class private functions
  ///////////////////////////////////////////////////////////////////////////////
  void        init();
  void        cleanup();
  int         next_non_empty_irrep(int n);
  ///////////////////////////////////////////////////////////////////////////////
  // Class data
  ///////////////////////////////////////////////////////////////////////////////
  // Type                           // Name
  CCIndex*                          ccindex_;
  int                               nelements_;
  int**                             element_irrep;
  short**                           tuples;
  size_t                            max_index_in_irrep_;
  int                               select;
protected:
  static int                        nirreps_;
};

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_ccindex_iterator_h
