#include <stdio.h>
#include "psio.h"

void psio_tocprint(ULI unit, FILE *output)
{
  psio_tocentry *this_entry;

  this_entry = psio_unit[unit].toc;

  fprintf(output, "\nTable of Contents for Unit %5u\n", unit);
  fprintf(output,
"----------------------------------------------------------------------------\n");
  fprintf(output,
"Key                                   Spage    Soffset      Epage    Eoffset\n");
  fprintf(output,
"----------------------------------------------------------------------------\n");

  while(this_entry != NULL) {
      fprintf(output, "%-32s %10u %10u %10u %10u\n",
	      this_entry->key,
	      this_entry->sadd.page, this_entry->sadd.offset,
	      this_entry->eadd.page, this_entry->eadd.offset);
      this_entry = this_entry->next;
    }
  fprintf(output, "\n");
  fflush(output);
}
