
/* $Log$
 * Revision 1.1  2000/02/04 22:53:19  evaleev
 * Initial revision
 *
/* Revision 2.5  1999/11/01 20:10:55  evaleev
/* Added explicit extern declarations of functions within the library.
/*
/* Revision 2.4  1997/08/25 21:49:50  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 2.3  1997/06/23  12:25:46  crawdad
 * Multiple changes to libciomr: Moved "param.h" to "iomrparam.h" to avoid
 *     conflicts with similarly named system file under linux.  Corrected type
 *    casting in rread(), rwrit(), sread(), and swrit() functions.  Corrected
 *     unclosed tmpdisks.dat file in sequential.c.  Corrected block_matrix() to
 *    avoid malloc'ing zero-length arrays.
 *
 * -Daniel
 *
 * Revision 2.2  1991/09/18  20:47:23  seidl
 * dec changes
 *
 * Revision 2.1  1991/06/15  18:29:14  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id$";


#include "iomrparam.h"
#include "includes.h"

extern int io_locate(FILE *, char[]);

#define MAX_SEGMENT 10

#define DEBUG 0
#define REGEX 1
#define NOFORTRAN 1

void print_segments();
void token_to_segments(), free_segments();
int matching_segments();

int
get_param(token,format,val)
char *token,*format;
#ifdef DEC
char *val;
#else
void *val;
#endif
{
  int i;
  char loc_token[MAX_STRING];
  char line[MAX_STRING];
  char first[MAX_STRING];
  char newformat[MAX_STRING];
  char *ptr;
  FILE *input;
  int ierr;
  FILE *fd;
  char *token_segments[MAX_SEGMENT];
  char *first_segments[MAX_SEGMENT];

#if DEBUG
  fprintf(stderr,"get_param: token: %s\n",token);
  (stderr,"get_param: format: %s\n",format);
#endif

  /* Initialize token_segments and first_segments to NULL's. */
  for (i=0; i<MAX_SEGMENT; i++) {
    token_segments[i] = NULL;
    first_segments[i] = NULL;
    }

  /* Build a string that locate can use. */
  strcpy(loc_token,"# ");
  strcat(loc_token,token);
  ptr = strchr(loc_token,':');
  strcpy(ptr," ##########");
  loc_token[10] = '\0';

  /* Open input.dat as the input file. */
  input = fopen("input.dat","r");
  if (input == NULL) {
    fprintf(stdout,"get_param: fopen: could not open input.dat\n");
    fprintf(stderr,"get_param: fopen: could not open input.dat\n");
    }
  ierr = io_locate(input,loc_token);
  if (ierr != 0) {
    fprintf(stdout,"get_param: io_locate: ierr=%d, token=%s\n", ierr,loc_token);
    fprintf(stderr,"get_param: io_locate: ierr=%d, token=%s\n", ierr,loc_token);
    return(-1);
    }

  ptr = strchr(token,':') + 1;

  /* Break the token into its segments. */
  token_to_segments(ptr,token_segments);

#if DEBUG
    fprintf(stderr,"get_param: ptr: %s\n",ptr);
#endif
  if (io_getline(input,line)!=0) return(-1);
  for (; line[0]!='#'; ) {
    /* Truncate line at 80 characters. */
    line[80] = '\0';
    sscanf(line,"%s",first);
    /* DEC compiler requires the following. */
    if (line[0] == '\0') first[0] = '\0';
    /* Break first into its segments. */
    token_to_segments(first,first_segments);
#if DEBUG
    print_segments(first_segments);
    fprintf(stderr,"get_param: line (1): %s\n",line);
    fprintf(stderr,"get_param: ptr: %s\n",ptr);
    fprintf(stderr,"get_param: first: %s\n",first);
#endif
    if (matching_segments(first_segments,token_segments)) {
      strcpy(newformat,"%s ");
      strcat(newformat,format);
      sscanf(line,newformat,first,val);
#if DEBUG
      fprintf(stderr,"get_param: line (2): %s\nget_param: data: ",line);
      fprintf(stderr,newformat,first,val);
      fprintf(stderr,"\n");
#endif
      fclose(input);
      free_segments(token_segments);
      free_segments(first_segments);
      return(0);
      }
    if (io_getline(input,line)!=0) return(-1);
    free_segments(first_segments);
    }
  fclose(input);
  free_segments(token_segments);
  return(-1);
  }

void
token_to_segments(token,seg)
char *token;
char *seg[MAX_SEGMENT];
{
  char *tloop;
  char *sloop;
  int  segn;
  int i;

  for (i=0; i<MAX_SEGMENT; i++) {
    seg[i] = NULL;
    }

  token[strlen(token)+1] = '\0';
  token[strlen(token)] = ':';

  sloop = token;
  segn = 0;
  for (tloop=token; *tloop!='\0'; tloop++) {
    if (segn >= MAX_SEGMENT) {
      fprintf(stderr,"token_to_segments: too many segments\n");
      exit(3);
      }
    else if (*tloop == ':') {
      seg[segn] = (char *) malloc(tloop - sloop + 1);
      strncpy(seg[segn],sloop,tloop - sloop);
      seg[segn][tloop - sloop] = '\0';
      segn++;
      sloop = tloop + 1;
      }
    }
  token[strlen(token)-1] = '\0';
  }

void
free_segments(seg)
char *seg[MAX_SEGMENT];
{
  int i;
  for (i=0; i<MAX_SEGMENT; i++) {
    if (!seg[i]) {
      free(seg[i]);
      seg[i] = NULL;
      }
    }
  }

int
matching_segments(seg1,seg2)
char *seg1[MAX_SEGMENT];
char *seg2[MAX_SEGMENT];
{
  int i;

  for (i=0; i<MAX_SEGMENT; i++) {
#if DEBUG
    fprintf(stderr,"matching_segments: %s %s\n",seg1[i],seg2[i]);
#endif
    if ((seg1[i] == NULL)&&(seg2[i] == NULL)) return(1);
    if (!matching_segment(seg1[i],seg2[i])) return(0);
    }
  return(1);
  }

int
matching_segment(seg1,seg2)
char *seg1, *seg2;
{
#if REGEX
  char recompch;
  int result;
#endif

  if ((seg1 == NULL)&&(seg2 == NULL)) return(1);
  if (seg1 == NULL) return(0);
  if (seg2 == NULL) return(0);
  if ((seg1[0] == '\0')&&(seg2[0] == '\0')) return(1);
  if (seg1[0] == '\0') return(0);
#if REGEX
  recompch = re_comp(seg1);
  if (recompch) {
    fprintf(stderr,"regex: failed, string = %s\n",seg1);
    fprintf(stderr,"regex: %s\n",recompch);
    exit(1);
    }
  result = re_exec(seg2);
  if (result == -1) {
    fprintf(stderr,"regex: internal error: %s %s\n",seg1,seg2);
    exit(1);
    }
  return(result);
#else
  if (!(   (!strcmp(seg1,"*"))
         ||(!strcmp(seg1,seg2)))) return(0);
  return(1);
#endif
  }

void
print_segments(seg)
char *seg[MAX_SEGMENT];
{
  int i;

  fprintf(stderr,"print_segments: ");

  for (i=0; seg[i] != NULL; i++) {
    fprintf(stderr,"%s ",seg[i]);
    }
  fprintf(stderr,"\n");
  }
