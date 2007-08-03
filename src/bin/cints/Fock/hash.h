/*! \file 
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/
/*-----------------------------------
  Declarations of htable_entry, etc.
 -----------------------------------*/
#include <psitypes.h>
typedef struct {
    PSI_INT_LEAST64 key;
    int si, sj, sk, sl;
    double q4ijkl, q4ikjl, q4ilkj;
} htable_entry;

typedef struct {
 htable_entry *table;
 int size;
} htable_t;

#define EMPTY_KEY -1

PSI_INT_LEAST64 compute_key(int si, int sj, int sk, int sl);
void init_htable(htable_t *htable, int nirreps);
void free_htable(htable_t *htable);
int put_entry(htable_t *htable, PSI_INT_LEAST64 key,
              int si, int sj, int sk, int sl,
	      double q4ijkl, double q4ikjl, double q4iljk);
int find_entry(htable_t *htable, PSI_INT_LEAST64 key);



