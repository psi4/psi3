/*!
  \file strncpy.c
  By Edward Valeev
  \ingroup (QT)
*/

#include <string.h>
#include <libqt/qt.h>

char*
psi::strncpy(char* dest, const char* source, size_t n) {
  if (n > 0) {
    ::strncpy(dest,source,n);
  }
  dest[n-1] = '\0';
  return dest;
}
