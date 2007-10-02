/*! \file gprgid.cc
  \ingroup (CINTS)
  \brief program id
*/
//namespace psi { namespace CINTS {
extern "C" {
  char *gprgid()
  {
    char *prgid = "cints";
    
    return(prgid);
  }
}
//}
//}
