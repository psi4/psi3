
/*!
  \file slaterd.c
  Edward Valeev, June 2002
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include "slaterd.h"


#define PSIO_INIT if (!psio_state()) { \
    psio_init(); \
    need_to_init_psio = 1; \
  }

#define PSIO_OPEN(u,n) if (!psio_open_check(u)) { \
    psio_open((u),n); \
    unit_opened = 0; \
  }

#define PSIO_CLOSE(u) if (!unit_opened) \
    psio_close((u),1);

#define PSIO_DONE if (need_to_init_psio) \
    psio_done();


/*! stringset_init()
 */
void stringset_init(StringSet *sset, int size, int nelec, int nfzc)
{
  sset->size = size;
  sset->nelec = nelec;
  sset->nfzc = nfzc;
  sset->strings = (String *) malloc(size*sizeof(String));
  memset(sset->strings,0,size*sizeof(String));
}

/*! stringset_delete()
 */
void stringset_delete(StringSet *sset)
{
  sset->size = 0;
  sset->nelec = 0;
  sset->nfzc = 0;
  if (sset->strings) free(sset->strings);
  sset->strings = NULL;
}

/*! stringset_add
 */
void stringset_add(StringSet *sset, int index, unsigned char *Occ)
{
  int i;
  int nact = sset->nelec - sset->nfzc;
  String *s;

  if (index < sset->size && index >= 0) {
    s = sset->strings + index;
  }
  s->index = index;
  s->occ = (short int*) malloc(nact*sizeof(short int));
  for(i=0;i<nact;i++)
    s->occ[i] = Occ[i];
}

/*! stringset_reindex
 */
void stringset_reindex(StringSet* sset, short int* mo_map)
{
  int s, mo;
  short int* occ;
  int nstrings = sset->size;
  int nact = sset->nelec - sset->nfzc;

  for(s=0; s<nstrings; s++) {
    occ = (sset->strings + s)->occ;
    for(mo=0; mo<nact; mo++)
      occ[mo] = mo_map[occ[mo]];
  }
}

void stringset_write(ULI unit, char *prefix, StringSet *sset)
{
  int i, size, nact;
  int need_to_init_psio = 0;
  int unit_opened = 1;
  char *size_key, *nelec_key, *nfzc_key, *strings_key;
  psio_address ptr;

PSIO_INIT
PSIO_OPEN(unit,PSIO_OPEN_OLD)

  size_key = (char *) malloc( strlen(prefix) + strlen(STRINGSET_KEY_SIZE) + 3);
  sprintf(size_key,":%s:%s",prefix,STRINGSET_KEY_SIZE);
  nelec_key = (char *) malloc( strlen(prefix) + strlen(STRINGSET_KEY_NELEC) + 3);
  sprintf(nelec_key,":%s:%s",prefix,STRINGSET_KEY_NELEC);
  nfzc_key = (char *) malloc( strlen(prefix) + strlen(STRINGSET_KEY_NFZC) + 3);
  sprintf(nfzc_key,":%s:%s",prefix,STRINGSET_KEY_NFZC);
  strings_key = (char *) malloc( strlen(prefix) + strlen(STRINGSET_KEY_STRINGS) + 3);
  sprintf(strings_key,":%s:%s",prefix,STRINGSET_KEY_STRINGS);

  psio_write_entry( unit, size_key, (char *)&sset->size, sizeof(int));
  psio_write_entry( unit, nelec_key, (char *)&sset->nelec, sizeof(int));
  psio_write_entry( unit, nfzc_key, (char *)&sset->nfzc, sizeof(int));
  ptr = PSIO_ZERO;
  size = sset->size;
  nact = sset->nelec - sset->nfzc;
  for(i=0; i<size; i++) {
    psio_write( unit, strings_key, (char *) &(sset->strings[i].index), sizeof(int), ptr, &ptr);
    psio_write( unit, strings_key, (char *) sset->strings[i].occ, nact*sizeof(short int), ptr, &ptr);
  }

PSIO_CLOSE(unit)
PSIO_DONE

  free(size_key);
  free(nelec_key);
  free(nfzc_key);
  free(strings_key);
}

void stringset_read(ULI unit, char *prefix, StringSet **stringset)
{
  int i, size, nelec, nfzc, nact;
  int need_to_init_psio = 0;
  int unit_opened = 1;
  char *size_key, *nelec_key, *nfzc_key, *strings_key;
  psio_address ptr;
  StringSet *sset = (StringSet *) malloc(sizeof(StringSet));

PSIO_INIT
PSIO_OPEN(unit,PSIO_OPEN_OLD)

  size_key = (char *) malloc( strlen(prefix) + strlen(STRINGSET_KEY_SIZE) + 3);
  sprintf(size_key,":%s:%s",prefix,STRINGSET_KEY_SIZE);
  nelec_key = (char *) malloc( strlen(prefix) + strlen(STRINGSET_KEY_NELEC) + 3);
  sprintf(nelec_key,":%s:%s",prefix,STRINGSET_KEY_NELEC);
  nfzc_key = (char *) malloc( strlen(prefix) + strlen(STRINGSET_KEY_NFZC) + 3);
  sprintf(nfzc_key,":%s:%s",prefix,STRINGSET_KEY_NFZC);
  strings_key = (char *) malloc( strlen(prefix) + strlen(STRINGSET_KEY_STRINGS) + 3);
  sprintf(strings_key,":%s:%s",prefix,STRINGSET_KEY_STRINGS);

  psio_read_entry( unit, size_key, (char *)&size, sizeof(int));
  psio_read_entry( unit, nelec_key, (char *)&nelec, sizeof(int));
  psio_read_entry( unit, nfzc_key, (char *)&nfzc, sizeof(int));

  stringset_init(sset, size, nelec, nfzc);

  nact = nelec - nfzc;
  ptr = PSIO_ZERO;
  for(i=0; i<size; i++) {
    psio_read( unit, strings_key, (char *) &(sset->strings[i].index), sizeof(int), ptr, &ptr);
    sset->strings[i].occ = (short int*) malloc(nact*sizeof(short int));
    psio_read( unit, strings_key, (char *) sset->strings[i].occ, nact*sizeof(short int), ptr, &ptr);
  }

PSIO_CLOSE(unit)
PSIO_DONE

  free(size_key);
  free(nelec_key);
  free(nfzc_key);
  free(strings_key);

  *stringset = sset;
}


/*! slaterdetset_init()
 */
void slaterdetset_init(SlaterDetSet *sdset, int size, StringSet *alphastrings, StringSet *betastrings)
{
  sdset->size = size;
  sdset->dets = (SlaterDet *) malloc(size*sizeof(SlaterDet));
  memset(sdset->dets,0,size*sizeof(SlaterDet));
  sdset->alphastrings = alphastrings;
  sdset->betastrings = betastrings;
}

/*! slaterdetset_delete()
 */
void slaterdetset_delete(SlaterDetSet *sdset)
{
  sdset->size = 0;
  if (sdset->dets) {
    free(sdset->dets);
    sdset->dets = NULL;
  }
  sdset->alphastrings = NULL;
  sdset->betastrings = NULL;
}

/*! slaterdetset_delete_full()
 */
void slaterdetset_delete_full(SlaterDetSet *sdset)
{
  sdset->size = 0;
  if (sdset->dets) {
    free(sdset->dets);
    sdset->dets = NULL;
  }
  if (sdset->alphastrings) {
    stringset_delete(sdset->alphastrings);
    sdset->alphastrings = NULL;
  }
  if (sdset->betastrings) {
    stringset_delete(sdset->betastrings);
    sdset->betastrings = NULL;
  }
}

/*! slaterdetset_add
 */
void slaterdetset_add(SlaterDetSet *sdset, int index, int alphastring, int betastring)
{
  SlaterDet *det;
  StringSet *alphaset = sdset->alphastrings;
  StringSet *betaset = sdset->betastrings;

  if (index < sdset->size && index >= 0) {
    det = sdset->dets + index;
  }
  det->index = index;
  if (alphastring < alphaset->size && alphastring >= 0)
    det->alphastring = alphastring;
  if (betastring < betaset->size && betastring >= 0)
    det->betastring = betastring;
}

void slaterdetset_write(ULI unit, char *prefix, SlaterDetSet *sdset)
{
  int i;
  int need_to_init_psio = 0;
  int unit_opened = 1;
  char *size_key, *set_key;
  char *alphaprefix, *betaprefix;
  psio_address ptr;

PSIO_INIT
PSIO_OPEN(unit,PSIO_OPEN_OLD)

  alphaprefix = (char *) malloc( strlen(prefix) + strlen(SDSET_KEY_ALPHASTRINGS) + 2);
  sprintf(alphaprefix,"%s:%s",prefix,SDSET_KEY_ALPHASTRINGS);
  betaprefix = (char *) malloc( strlen(prefix) + strlen(SDSET_KEY_BETASTRINGS) + 2);
  sprintf(betaprefix,"%s:%s",prefix,SDSET_KEY_BETASTRINGS);

  stringset_write( unit, alphaprefix, sdset->alphastrings);
  stringset_write( unit, betaprefix, sdset->betastrings);
  
  free(alphaprefix);
  free(betaprefix);

  size_key = (char *) malloc( strlen(prefix) + strlen(SDSET_KEY_SIZE) + 3);
  sprintf(size_key,":%s:%s",prefix,SDSET_KEY_SIZE);
  set_key = (char *) malloc( strlen(prefix) + strlen(SDSET_KEY_DETERMINANTS) + 3);
  sprintf(set_key,":%s:%s",prefix,SDSET_KEY_DETERMINANTS);

  psio_write_entry( unit, size_key, (char *)&sdset->size, sizeof(int));
  psio_write_entry( unit, set_key, (char *)sdset->dets, sdset->size*sizeof(SlaterDet));

PSIO_CLOSE(unit)
PSIO_DONE

  free(size_key);
  free(set_key);
}

void slaterdetset_read(ULI unit, char *prefix, SlaterDetSet **slaterdetset)
{
  int i, size;
  int need_to_init_psio = 0;
  int unit_opened = 1;
  char *size_key, *set_key;
  char *alphaprefix, *betaprefix;
  psio_address ptr;
  StringSet *alphastrings, *betastrings;
  SlaterDetSet *sdset = (SlaterDetSet *) malloc(sizeof(SlaterDetSet));

PSIO_INIT
PSIO_OPEN(unit,PSIO_OPEN_OLD)

  alphaprefix = (char *) malloc( strlen(prefix) + strlen(SDSET_KEY_ALPHASTRINGS) + 2);
  sprintf(alphaprefix,"%s:%s",prefix,SDSET_KEY_ALPHASTRINGS);
  betaprefix = (char *) malloc( strlen(prefix) + strlen(SDSET_KEY_BETASTRINGS) + 2);
  sprintf(betaprefix,"%s:%s",prefix,SDSET_KEY_BETASTRINGS);

  stringset_read( unit, alphaprefix, &alphastrings);
  stringset_read( unit, betaprefix, &betastrings);
  
  free(alphaprefix);
  free(betaprefix);

  size_key = (char *) malloc( strlen(prefix) + strlen(SDSET_KEY_SIZE) + 3);
  sprintf(size_key,":%s:%s",prefix,SDSET_KEY_SIZE);
  set_key = (char *) malloc( strlen(prefix) + strlen(SDSET_KEY_DETERMINANTS) + 3);
  sprintf(set_key,":%s:%s",prefix,SDSET_KEY_DETERMINANTS);

  psio_read_entry( unit, size_key, (char *)&size, sizeof(int));
  slaterdetset_init(sdset,size,alphastrings,betastrings);
  psio_read_entry( unit, set_key, (char *)sdset->dets, sdset->size*sizeof(SlaterDet));

PSIO_CLOSE(unit)
PSIO_DONE

  free(size_key);
  free(set_key);

  *slaterdetset = sdset;
}


/*! slaterdetvector_init()
 */
void slaterdetvector_init(SlaterDetVector *sdvector, SlaterDetSet *sdset)
{
  sdvector->size = sdset->size;
  sdvector->sdset = sdset;
  sdvector->coeffs = init_array(sdvector->size);
}

/*! slaterdetvector_delete()
 */
void slaterdetvector_delete(SlaterDetVector *sdvector)
{
  sdvector->size = 0;
  sdvector->sdset = NULL;
  if (sdvector->coeffs) {
    free(sdvector->coeffs);
    sdvector->coeffs = NULL;
  }
}

/*! slaterdetvector_delete_full()
 */
void slaterdetvector_delete_full(SlaterDetVector *sdvector)
{
  sdvector->size = 0;
  if (sdvector->sdset) {
    slaterdetset_delete_full(sdvector->sdset);
    sdvector->sdset = NULL;
  }
  if (sdvector->coeffs) {
    free(sdvector->coeffs);
    sdvector->coeffs = NULL;
  }
}

/*! slaterdetvector_add
 */
void slaterdetvector_add(SlaterDetVector *sdvector, int index, double coeff)
{
  if (index < sdvector->size && index >= 0) {
    sdvector->coeffs[index] = coeff;
  }
}

/*! slaterdetvector_set
 */
void slaterdetvector_set(SlaterDetVector *sdvector, double *coeffs)
{
  int i;
  const int size = sdvector->size;
  double *v = sdvector->coeffs;
  if (v) {
    for(i=0; i<size; i++)
      v[i] = coeffs[i];
  }
}


void slaterdetvector_write(ULI unit, char *prefix, SlaterDetVector *vector)
{
  int i;
  int need_to_init_psio = 0;
  int unit_opened = 1;
  char *vector_key;
  psio_address ptr;

PSIO_INIT
PSIO_OPEN(unit,PSIO_OPEN_OLD)

  slaterdetset_write( unit, prefix, vector->sdset);

  vector_key = (char *) malloc( strlen(prefix) + strlen(SDVECTOR_KEY_VECTOR) + 3);
  sprintf(vector_key,":%s:%s",prefix,SDVECTOR_KEY_VECTOR);

  psio_write_entry( unit, vector_key, (char *)vector->coeffs, vector->size*sizeof(double));

PSIO_CLOSE(unit)
PSIO_DONE

  free(vector_key);
}


void slaterdetvector_read(ULI unit, char *prefix, SlaterDetVector **sdvector)
{
  int i;
  int need_to_init_psio = 0;
  int unit_opened = 1;
  char *vector_key;
  psio_address ptr;
  SlaterDetSet *sdset;
  SlaterDetVector *vector = (SlaterDetVector *) malloc(sizeof(SlaterDetVector));

PSIO_INIT
PSIO_OPEN(unit,PSIO_OPEN_OLD)

  slaterdetset_read( unit, prefix, &sdset);
  slaterdetvector_init(vector,sdset);

  vector_key = (char *) malloc( strlen(prefix) + strlen(SDVECTOR_KEY_VECTOR) + 3);
  sprintf(vector_key,":%s:%s",prefix,SDVECTOR_KEY_VECTOR);

  psio_read_entry( unit, vector_key, (char *)vector->coeffs, vector->size*sizeof(double));

PSIO_CLOSE(unit)
PSIO_DONE
  
  free(vector_key);

  *sdvector = vector;
}


