/*! \file
  \ingroup CINTS
  \brief program id
*/
extern "C" {
  const char *gprgid()
  {
    const char *prgid = "cints";
    
    return(prgid);
  }
}
