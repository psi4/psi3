/*!
 ** \file
 ** \ingroup PSIO
 */

#include <cstdlib>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

#ifdef PSIO_STATS
#include <<ctime>>
#endif

using namespace psi;

PSIO::~PSIO() {
#ifdef PSIO_STATS
  int i;
  ULI total_read=0, total_write=0;
  FILE *io_out;
  time_t my_time;
  my_time = time(NULL);
  io_out = fopen("io.dat", "a+");
  fprintf(io_out, "\nLIBPSIO Read/Write Statistics\n\n");
  fprintf(io_out, "Run at: %s\n", ctime(&my_time));
  fprintf(io_out, "Unit      Read(kB)    Write(kB)\n");
  fprintf(io_out, "-------------------------------\n");
  for(i=0; i < psio_unit.size(); i++) {
    total_read += psio_unit[i].readlen;
    total_write += psio_unit[i].writlen;

    if(psio_unit[i].readlen || psio_unit[i].writlen)
    fprintf(io_out, "%3d   %10.1f   %10.1f\n",i,
        ((double) psio_unit[i].readlen)/((double) 1024),
        ((double) psio_unit[i].writlen)/((double) 1024));
  }
  fprintf(io_out, "-------------------------------\n");
  fprintf(io_out, "Total %10.1f   %10.1f\n",
      ((double) total_read)/((double) 1024),
      ((double) total_write)/((double) 1024));
  fclose(io_out);
#endif
  
  state_ = 0;
}

extern "C" {
  /*!
   ** \ingroup PSIO
   **
   ** PSIO_DONE(): Frees global data used by the I/O routines.
   **
   ** No arguments.
   */

  int psio_done(void) {
    delete _default_psio_lib_;
  }

}
