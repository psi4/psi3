#include <stdio.h>
#include <stdlib.h>
#include <qt.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

void dpd_file4_cache_init(void)
{
  dpd_default->file4_cache = NULL;
  dpd_default->file4_cache_most_recent = 0;
  dpd_default->file4_cache_least_recent = 1;
}

void dpd_file4_cache_close(void)
{
  struct dpd_file4_cache_entry *this_entry, *next_entry;
  dpdfile4 Outfile;

  this_entry = dpd_default->file4_cache;
  
  while(this_entry != NULL) {

      /* Clean out each file4_cache entry */
      dpd_file4_init(&Outfile, this_entry->filenum, this_entry->irrep,
		    this_entry->pqnum, this_entry->rsnum, this_entry->label);

      next_entry = this_entry->next;
      
      dpd_file4_cache_del(&Outfile);
      dpd_file4_close(&Outfile);

      this_entry = next_entry;
   }
}

struct dpd_file4_cache_entry 
*dpd_file4_cache_scan(int filenum, int irrep, int pqnum, int rsnum, char *label)
{
  struct dpd_file4_cache_entry *this_entry;

#ifdef DPD_TIMER
  timer_on("file4_cache");
#endif

  this_entry = dpd_default->file4_cache;

  while(this_entry != NULL) {
      if(this_entry->filenum == filenum      &&
	 this_entry->irrep == irrep          &&
         this_entry->pqnum == pqnum          &&
         this_entry->rsnum == rsnum          &&
         !strcmp(this_entry->label,label)) {
#ifdef DPD_TIMER
              timer_off("file4_cache");
#endif
	      /* increment the access timers */
	      dpd_default->file4_cache_most_recent++;
	      this_entry->access = dpd_default->file4_cache_most_recent;
	      
              return(this_entry);
           }
       
      this_entry = this_entry->next;
    }

#ifdef DPD_TIMER
  timer_off("file4_cache");
#endif
  return(this_entry);
}

struct dpd_file4_cache_entry *dpd_file4_cache_last(void)
{
  struct dpd_file4_cache_entry *this_entry;

  this_entry = dpd_default->file4_cache;

  while(this_entry !=NULL) {
      if(this_entry->next == NULL) return(this_entry);
      this_entry = this_entry->next;
    }

  return(NULL);
}

int dpd_file4_cache_add(dpdfile4 *File)
{
  int h;
  struct dpd_file4_cache_entry *this_entry;

  this_entry = dpd_file4_cache_scan(File->filenum, File->my_irrep,
				   File->params->pqnum, File->params->rsnum,
                                   File->label);

  if((this_entry != NULL && !(File->incore)) ||
     (this_entry == NULL && (File->incore))) {
      /* Either the file4 appears in the cache but incore is not set,
	 or incore is set and the file4 isn't in the cache */
      dpd_error("File4 cache add error!", stderr);
    }
  else if(this_entry != NULL && File->incore) 
      return 0; /* We already have this one in cache */
  else if(this_entry == NULL && !(File->incore)) { /* New cache entry */
      this_entry = (struct dpd_file4_cache_entry *) 
                       malloc(sizeof(struct dpd_file4_cache_entry));
      this_entry->filenum = File->filenum;
      this_entry->irrep = File->my_irrep;
      this_entry->pqnum = File->params->pqnum;
      this_entry->rsnum = File->params->rsnum;
      strcpy(this_entry->label,File->label);
      this_entry->next = NULL;
      this_entry->last = dpd_file4_cache_last();
      
      if(this_entry->last != NULL) this_entry->last->next = this_entry;
      else dpd_default->file4_cache = this_entry;

      /* increment the access timers */
      dpd_default->file4_cache_most_recent++;
      this_entry->access = dpd_default->file4_cache_most_recent;
	      
      /* Read all data into core */
      this_entry->size = 0;
      for(h=0; h < File->params->nirreps; h++) {
	  this_entry->size += 
	      File->params->rowtot[h] * File->params->coltot[h];
	  dpd_file4_mat_irrep_init(File, h);
	  dpd_file4_mat_irrep_rd(File, h);
	}

      this_entry->matrix = File->matrix;

      File->incore = 1;

      return 0;
  }

  return 0;
}

int dpd_file4_cache_del(dpdfile4 *File)
{
  int h;
  struct dpd_file4_cache_entry *this_entry, *next_entry, *last_entry;

  this_entry = dpd_file4_cache_scan(File->filenum, File->my_irrep,
				    File->params->pqnum, File->params->rsnum,
				    File->label);

  if((this_entry == NULL && File->incore) ||
     (this_entry != NULL && !(File->incore)) ||
      (this_entry == NULL && !(File->incore))) {
      dpd_error("File4 cache delete error!", stderr);
    }
  else {
      File->incore = 0;

      /* Write all the data to disk and free the memory */
      for(h=0; h < File->params->nirreps; h++) {
	  dpd_file4_mat_irrep_wrt(File, h);
	  dpd_file4_mat_irrep_close(File, h);
	}

      next_entry = this_entry->next;
      last_entry = this_entry->last;

      /* Are we deleting the top of the tree? */
      if(this_entry == dpd_default->file4_cache) 
	  dpd_default->file4_cache = next_entry;

      free(this_entry);

      /* Reassign pointers for adjacent entries in the list */
      if(next_entry != NULL) next_entry->last = last_entry;
      if(last_entry != NULL) last_entry->next = next_entry;

    }

  return 0;
}

void dpd_file4_cache_print(FILE *outfile)
{
  int total_size=0;
  struct dpd_file4_cache_entry *this_entry;

  this_entry = dpd_default->file4_cache;

  fprintf(outfile, "\n\tDPD File4 Cache Listing:\n\n");
  fprintf(outfile,
    "Cache Label                     File symm  pq  rs  access  size(kB)\n");
  fprintf(outfile,
    "-------------------------------------------------------------------\n");
  while(this_entry != NULL) {
      fprintf(outfile,
      "%-32s %3d    %1d  %2d  %2d  %6d  %8.1f\n",
	      this_entry->label, this_entry->filenum, this_entry->irrep,
	      this_entry->pqnum, this_entry->rsnum, this_entry->access,
	      (this_entry->size)*sizeof(double)/1e3);
      total_size += this_entry->size;
      this_entry = this_entry->next;
    }
  fprintf(outfile,
    "-----------------------------------------------------------\n");
  fprintf(outfile, "Total cached: %8.1f kB; MRU = %6d; LRU = %6d\n",
	  (total_size*sizeof(double))/1e3, dpd_default->file4_cache_most_recent,
	  dpd_default->file4_cache_least_recent);
}

struct dpd_file4_cache_entry *dpd_file4_cache_find_lru(void)
{
  struct dpd_file4_cache_entry *this_entry;

  this_entry = dpd_default->file4_cache;

  if(this_entry == NULL) return(NULL);

  while(dpd_default->file4_cache_least_recent <=
	dpd_default->file4_cache_most_recent) {
      while(this_entry !=NULL) {
	  if(this_entry->access == dpd_default->file4_cache_least_recent)
	      return(this_entry);
	  this_entry = this_entry->next;
	}
      dpd_default->file4_cache_least_recent++;
      this_entry = dpd_default->file4_cache;
    }

  dpd_error("Error locating file4_cache LRU!", stderr);
}

int dpd_file4_cache_del_lru(void)
{
  dpdfile4 File;
  struct dpd_file4_cache_entry *this_entry;

  this_entry = dpd_file4_cache_find_lru();

  if(this_entry == NULL) return 1; /* there is no cache */
  else { /* we found the LRU so delete it */
      printf("Deleteing LRU: %-32s %3d %2d %2d %6d\n", this_entry->label,
	     this_entry->filenum, this_entry->pqnum, this_entry->rsnum,
	     this_entry->access);
      
      dpd_file4_init(&File, this_entry->filenum, this_entry->irrep,
		     this_entry->pqnum, this_entry->rsnum, this_entry->label);

      dpd_file4_cache_del(&File);
      dpd_file4_close(&File);

      return 0;
    }
}
