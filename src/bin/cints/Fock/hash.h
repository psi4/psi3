/*-----------------------------------
  Declarations of htable_entry, etc.
 -----------------------------------*/
typedef struct {
    long int key;
    int si, sj, sk, sl;
    double q4ijkl, q4ikjl, q4ilkj;
} htable_entry;

typedef struct {
 htable_entry *table;
 int size;
} htable_t;

#define EMPTY_KEY -1

long int compute_key(int si, int sj, int sk, int sl);
void init_htable(htable_t *htable, int nirreps);
void free_htable(htable_t *htable);
int put_entry(htable_t *htable, long int key, int si, int sj, int sk, int sl,
	      double q4ijkl, double q4ikjl, double q4iljk);
int find_entry(htable_t *htable, long int key);



