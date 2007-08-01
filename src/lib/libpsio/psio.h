#ifndef PSIO_H
#define PSIO_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#define PSIO_OPEN_NEW 0
#define PSIO_OPEN_OLD 1

#define PSIO_KEYLEN 80
#define PSIO_MAXVOL 8
#define PSIO_MAXUNIT 300
#define PSIO_PAGELEN 65536

typedef unsigned long int ULI;  /* For convenience */

typedef struct {
    ULI page;   /* First page of entry */
    ULI offset; /* Starting byte offset on fpage */
} psio_address;

struct psio_entry {
    char key[PSIO_KEYLEN];
    psio_address sadd;
    psio_address eadd;
    struct psio_entry *next;
    struct psio_entry *last;
};

typedef struct psio_entry psio_tocentry;

typedef struct {
    char *path;
    int stream;
} psio_vol;

typedef struct {
    ULI numvols;
    psio_vol vol[PSIO_MAXVOL];
    ULI toclen;
    psio_tocentry *toc;
} psio_ud;

/**
  default libpsio is globally visible
*/
typedef struct {
  psio_ud *psio_unit;
#ifdef PSIO_STATS
  ULI *psio_readlen;
  ULI *psio_writlen;
#endif
  /** Library state variable */
  int state;
} psio_lib;
/** the default library */
extern psio_lib* _default_psio_lib_;
/** Upon catastrophic failure, the library will exit() with this code. The default is 1, but can be overridden. */
extern int _psio_error_exit_code_;
/** A convenient address initialization struct */
extern psio_address PSIO_ZERO;

#define PSIO_ERROR_INIT       1
#define PSIO_ERROR_DONE       2
#define PSIO_ERROR_MAXVOL     3
#define PSIO_ERROR_NOVOLPATH  4
#define PSIO_ERROR_OPEN       5
#define PSIO_ERROR_REOPEN     6
#define PSIO_ERROR_CLOSE      7
#define PSIO_ERROR_RECLOSE    8
#define PSIO_ERROR_OSTAT      9
#define PSIO_ERROR_LSEEK     10
#define PSIO_ERROR_READ      11
#define PSIO_ERROR_WRITE     12
#define PSIO_ERROR_NOTOCENT  13
#define PSIO_ERROR_TOCENTSZ  14
#define PSIO_ERROR_KEYLEN    15
#define PSIO_ERROR_BLKSIZ    16
#define PSIO_ERROR_BLKSTART  17
#define PSIO_ERROR_BLKEND    18
#define PSIO_ERROR_IDENTVOLPATH 19
#define PSIO_ERROR_MAXUNIT   20

int psio_init(void);
int psio_init_ipfree(void);
int psio_state(void);
int psio_done(void);
void psio_error(unsigned int unit, unsigned int errval);
int psio_open(unsigned int unit, int status);
int psio_close(unsigned int unit, int keep);

unsigned int psio_get_numvols(unsigned int unit);
unsigned int psio_get_numvols_default(void);
int psio_get_volpath(unsigned int unit, unsigned int volume, char **path);
int psio_get_volpath_default(unsigned int volume, char **path);
int psio_get_filename(unsigned int unit, char **name);
int psio_get_filename_default(char **name);
psio_address psio_get_address(psio_address start, ULI shift);
psio_address psio_get_global_address(psio_address entry_start,
                                     psio_address rel_address);
int psio_volseek(psio_vol *vol, ULI page, ULI offset, ULI numvols);
ULI psio_get_length(psio_address sadd, psio_address eadd);
psio_address psio_get_entry_end(unsigned int unit, char *key);

int psio_tocwrite(unsigned int unit);
int psio_tocread(unsigned int unit);
void psio_tocprint(unsigned int unit, FILE *output);
psio_tocentry *psio_tocscan(unsigned int unit, char *key);
psio_tocentry *psio_toclast(unsigned int unit);
unsigned int psio_toclen(unsigned int unit);
int psio_tocdel(unsigned int unit, char *key);
int psio_tocclean(unsigned int unit, char *key);
void psio_tocrename(unsigned int unit, char *key, char *newkey);

int psio_write(unsigned int unit, char *key, char *buffer, ULI size,
	       psio_address sadd, psio_address *eadd);
int psio_read(unsigned int unit, char *key, char *buffer, ULI size,
	      psio_address sadd, psio_address *eadd);
int psio_write_entry(unsigned int unit, char *key, char *buffer, ULI size);
int psio_read_entry(unsigned int unit, char *key, char *buffer, ULI size);
int psio_write_block(unsigned int unit, char *key, char *buffer, ULI blksiz,
		     ULI start_blk, ULI end_blk);
int psio_read_block(unsigned int unit, char *key, char *buffer, ULI blksiz,
		    ULI start_blk, ULI end_blk);
int psio_rw(unsigned int unit, char *buffer, psio_address address, ULI size, int wrt);

int psio_open_check(unsigned int unit);

ULI psio_rd_toclen(unsigned int unit);
void psio_wt_toclen(unsigned int unit, ULI toclen);

int psio_set_filescfg_kwd(const char* kwdgrp, const char* kwd, int unit, const char* kwdval);
const char* psio_get_filescfg_kwd(const char* kwdgrp, const char* kwd, int unit);

#ifdef __cplusplus
}
#endif

#endif    /* #ifndef PSIO_H */
