#ifndef PSIO_H
#define PSIO_H

#include <stdio.h>
#include "psio.gbl"

/* A convenient address initialization struct */
extern psio_address PSIO_ZERO;

#ifdef PSIO_STATS
extern ULI *psio_readlen;
extern ULI *psio_writlen;
#endif

int psio_init(void);
int psio_state(void);
int psio_done(void);
void psio_error(ULI unit, ULI errval);
int psio_open(ULI unit, int status);
int psio_close(ULI unit, int keep);

ULI psio_get_numvols(ULI unit);
ULI psio_get_numvols_default(void);
int psio_get_volpath(ULI unit, ULI volume, char *path);
int psio_get_volpath_default(ULI volume, char *path);
int psio_get_filename(ULI unit, char *name);
int psio_get_filename_default(char *name);
psio_address psio_get_address(psio_address start, ULI shift);
psio_address psio_get_global_address(psio_address entry_start,
                                     psio_address rel_address);
int psio_volseek(psio_vol *vol, ULI page, ULI offset, ULI numvols);
ULI psio_get_length(psio_address sadd, psio_address eadd);
psio_address psio_get_entry_end(ULI unit, char *key);

int psio_tocwrite(ULI unit);
int psio_tocread(ULI unit);
void psio_tocprint(ULI unit, FILE *output);
psio_tocentry *psio_tocscan(ULI unit, char *key);
psio_tocentry *psio_toclast(ULI unit);
ULI psio_toclen(ULI unit);
int psio_tocdel(ULI unit, char *key);
int psio_tocclean(ULI unit, char *key);
void psio_tocrename(ULI unit, char *key, char *newkey);

int psio_write(ULI unit, char *key, char *buffer, ULI size,
	       psio_address sadd, psio_address *eadd);
int psio_read(ULI unit, char *key, char *buffer, ULI size,
	      psio_address sadd, psio_address *eadd);
int psio_write_entry(ULI unit, char *key, char *buffer, ULI size);
int psio_read_entry(ULI unit, char *key, char *buffer, ULI size);
int psio_write_block(ULI unit, char *key, char *buffer, ULI blksiz,
		     ULI start_blk, ULI end_blk);
int psio_read_block(ULI unit, char *key, char *buffer, ULI blksiz,
		    ULI start_blk, ULI end_blk);
int psio_rw(ULI unit, char *buffer, psio_address address, ULI size, int wrt);

int psio_open_check(ULI unit);


#endif    /* #ifndef PSIO_H */
