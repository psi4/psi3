
#include <libpsio/psio.h>

#include <string>
#include <map>
#include <sstream>
#include <algorithm>
#include <cctype>

using namespace std;

namespace {
  typedef std::map<std::string,std::string> KWDMap;
  /// keeps keywords dealing with files configuration
  KWDMap files_keywords;

  std::string fullkwd(const char* kwdgrp, const char* kwd, int unit) {
    std::string unitname;
    if (unit < 0)
      unitname = "DEFAULT";
    else {
      std::ostringstream oss;
      oss << "FILE" << unit;
      unitname = oss.str();
    }
    const std::string sep(":");
 
    std::string fkwd = sep + kwdgrp + sep + "FILES" + sep + unitname + sep + kwd;
    // convert to upper case
    std::transform(fkwd.begin(), fkwd.end(), fkwd.begin(), 
		   static_cast<int(*)(int)>(toupper));
    return fkwd;
  }
}

extern "C" {

  /*!
    psio_set_filescfg_kwd(): set keyword kwd describing some aspect of configuration of PSIO file unit
    to value kwdval. kwdgrp specifies the keyword group (useful values are: "DEFAULT", "PSI", and the name of
    the current executable). If unit is set to -1, this keyword will set the default for all units (this keyword
    can be further overridden for some units). To specify a keyword that works for a specific unit, set unit to the
    appropriate number between 0 to PSIO_MAXUNIT.

    libpsio understands the following keywords: "name" (specifies the prefix for the filename,
    i.e. if name is set to "psi" then unit 35 will be named "psi.35"), "nvolume" (number of files over which
    to stripe this unit, cannot be greater than PSIO_MAXVOL), "volumeX", where X is a positive integer less than or equal to
    the value of "nvolume".

    \ingroup(PSIO)
  */
  int psio_set_filescfg_kwd(const char* kwdgrp, const char* kwd, int unit, const char* kwdval)
  {
    std::string fkwd = fullkwd(kwdgrp,kwd,unit);
    files_keywords[fkwd] = kwdval;
    return 0;
  }

  const char* psio_get_filescfg_kwd(const char* kwdgrp, const char* kwd, int unit)
  {
    static char* nullstr = 0;

    const std::string fkwd = fullkwd(kwdgrp,kwd,unit);
    KWDMap::const_iterator kwd_loc = files_keywords.find(fkwd);
    if (kwd_loc != files_keywords.end())
      return kwd_loc->second.c_str();
    else
      return nullstr;
  }

}
