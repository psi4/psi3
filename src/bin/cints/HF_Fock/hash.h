/*-----------------------------------
  Declarations of htable_entry, etc.
 -----------------------------------*/
typedef struct {
    int key;
    int si, sj, sk, sl;
    double q4ijkl, q4ikjl, q4ilkj;
} htable_entry;

int compute_key(int si, int sj, int sk, int sl);
void init_htable(int nirreps);
void free_htable();
int put_entry(int key, int si, int sj, int sk, int sl, double q4ijkl, double q4ikjl, double q4iljk);
int find_entry(int key);

#define EMPTY_KEY -1

extern htable_entry *htable;
