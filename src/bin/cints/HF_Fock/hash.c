#include<stdio.h>
#include<stdlib.h>
#include<libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include"hash.h"

#define hashing_function(a) (a)%htable_size   /* division method */

static int htable_size;
htable_entry *htable;

/*-------------------------
  Initialize hashing table
 -------------------------*/
void init_htable(int nirreps)
{
  int i;
  
  switch(nirreps) {
  case 1:
      fprintf(stderr," init_htable: nirreps = 1\n");
      exit(1);
  case 2:
      htable_size = 97;
      break;
  case 4:
      htable_size = 1543;
      break;
  case 8:
  default:
      htable_size = 24593;
      break;
  }

  if (htable == NULL) {
    htable = (htable_entry *) malloc(htable_size*sizeof(htable_entry));
    for(i=0;i<htable_size;i++)
      htable[i].key = EMPTY_KEY;
  }

  return;
}

      
/*-------------------------
  Deallocate hashing table
 -------------------------*/
void free_htable()
{
  free(htable);
  htable = NULL;
  return;
}


/*------------------------------------------
  Compute the key from four (shell) indices
 ------------------------------------------*/
int compute_key(int si, int sj, int sk, int sl)
{
  int ij, kl, ijkl;
  int dum;

/*  ij = (usi_eq_usj) ? INDEX(si,sj) : si*num_shells + sj;
  kl = (usk_eq_usl) ? INDEX(sk,sl) : sk*num_shells + sl;
  ijkl = (usij_eq_uskl) ? INDEX(ij,kl) : ij*num_shells*num_shells + kl;*/

  return INDEX(INDEX(si,sj),INDEX(sk,sl));
}


/*----------------------------------------------------
  int put_entry() :  Put an entry into the hashing table.
                     Return the location if it's a new entry,
		     or -1 if it has been in the table.

  One must remember that given set of si, sj, sk, sl and
  the set found in the table may be a permutation of each
  other which belong to Q4 - group which consists of all
  permutations under which an ERI is unchanged
  (Q4 is a subgroup of T4 of order 8).
  Another complication
  is that if si == sj or sk == sl - only htable[].q4ikjl
  is relevant, yet due to possibility of permutations
  it is safer to just increment it by THE SUM of q4ikjl and
  q4ilkj (since one of them is zero);

  Let P be the permutation which transforms the given {si,sj,sk,sl}
  into htable[].{si,sj,sk,sl}. P must be a product of the following
  (commuting) operations P0 (trivial transposition), P12 (swaps
  si and sj), P34 (swaps sk and sl), and P12,34 = P13 P24 (swaps bra
  and ket). Let K be a transposition that describes the contribution
  of a given {si,sj,sk,sl} to the current P(si,sj|sk,sl). K must be
  either P0 (q4ijkl is non-zero), P23 (q4ikjl non-zero) or P24
  (q4ilkj non-zero). The following set of rules may be easily derived
  using multiplication table of T4:

  P0  {P} equiv(Q4) P0
  P23 P12 equiv(Q4) P24
  P24 P12 equiv(Q4) P23
  P23 P34 equiv(Q4) P24
  P24 P34 equiv(Q4) P23
  P23 P12,34 equiv(Q4) P23
  P24 P12,34 equiv(Q4) P24

  equiv(Q4) means equivalent up to any permutation from Q4.

  This means that P is either P12 or P34 we need to swap q4ikjl and q4ilkj.

 ----------------------------------------------------*/
int put_entry(int key, int si, int sj, int sk, int sl, double q4ijkl, double q4ikjl, double q4ilkj)
{
  int curr_ptr;
  int hvalue = hashing_function(key);
  int return_code = -1;
  int P_includes_P12, P_includes_P34;
  
  curr_ptr = hvalue;
  while(htable[curr_ptr].key != key && htable[curr_ptr].key != EMPTY_KEY) {
    curr_ptr++;
    if (htable_size == curr_ptr)
      curr_ptr = 0;
  }
  if (htable[curr_ptr].key == EMPTY_KEY) {
    return_code = curr_ptr;
    htable[curr_ptr].key = key;
    htable[curr_ptr].si = si;
    htable[curr_ptr].sj = sj;
    htable[curr_ptr].sk = sk;
    htable[curr_ptr].sl = sl;
    htable[curr_ptr].q4ijkl = q4ijkl;
    if (si == sj || sk == sl) {
      htable[curr_ptr].q4ikjl = q4ikjl + q4ilkj;
      htable[curr_ptr].q4ilkj = 0;
    }
    else {
      htable[curr_ptr].q4ikjl = q4ikjl;
      htable[curr_ptr].q4ilkj = q4ilkj;
    }
  }
  else {
    htable[curr_ptr].q4ijkl += q4ijkl;
    if (si == sj || sk == sl) {
      htable[curr_ptr].q4ikjl += q4ikjl + q4ilkj;
    }
    else {
      P_includes_P34 = ( (htable[curr_ptr].sk == sl && htable[curr_ptr].sl == sk) ||
			 (htable[curr_ptr].sk == sj && htable[curr_ptr].sl == si) );
      P_includes_P12 = ( (htable[curr_ptr].si == sj && htable[curr_ptr].sj == si) ||
			 (htable[curr_ptr].si == sl && htable[curr_ptr].sj == sk) );
      if (P_includes_P12 ^ P_includes_P34) {
	htable[curr_ptr].q4ikjl += q4ilkj;
	htable[curr_ptr].q4ilkj += q4ikjl;
      }
      else {
	htable[curr_ptr].q4ikjl += q4ikjl;
	htable[curr_ptr].q4ilkj += q4ilkj;
      }
    }
  }

  return return_code;
}

/*----------------------------------------------------
  Put an entry into the table. Return the location if
  it's a new entry, or -1 if it has been in the table
 ----------------------------------------------------*/
int find_entry(int key)
{
  int curr_ptr;
  int hvalue = hashing_function(key);
  int return_code = -1;
  
  curr_ptr = hvalue;
  while(htable[curr_ptr].key != key && htable[curr_ptr].key != EMPTY_KEY) {
    curr_ptr++;
    if (htable_size == curr_ptr)
      curr_ptr = 0;
  }
  if (htable[curr_ptr].key == EMPTY_KEY)
    return -1;
  else
    return curr_ptr;

}
