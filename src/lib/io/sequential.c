
#include "includes.h"
#include "param.h"
#include "types.h"

#if defined(SUN)
#include <unistd.h>
#endif

#include <errno.h>

#if defined(DEC) || (defined(SGI) && !defined(__GNUG__))
extern char *getenv();
#endif

#ifndef FILENAME_MAX
#define FILENAME_MAX 512
#endif

/* CDS 9/93 */
#define TMP_VOL_MAX 10          /* max number of tmp drives */

#if defined(DEC)                
char *strstr() ;                /* this function not in old DEC compiler! */
#endif
/* CDS 9/93 */

#define DEBUG 0

static char open_name[] = "sequential_ioopen";

sequential_t *
sequential_ioopen(baseparam,unit)
char *baseparam;
int unit;
{
  int osi,bool;
  int i;
  sequential_t *ud;
  char param[MAX_STRING];
  char volid[MAX_STRING];
  char name[MAX_STRING];
  char volpath[MAX_STRING];
  char path[MAX_STRING];
  char *mode;
  struct stat statjunk;
  char name_format[3];
  /* This file contains a list of files to be deleted
   * when we are using batch. */
  FILE *cleanfilep;
  /* If we are running in batch then this is set to the $(MBATCH)
   * environment variable to indicate the job_id. */
  int job_id;
  char *mbatchc;
  char cleanfile[FILENAME_MAX];

/* CDS 9/93 New variables for automatic tmp drive determination */
  char endofvolpath[MAX_STRING] ; /* last part of volume path string */
  int got_first_volpath ;         /* flag for getting first volume path */
  int use_default_vols=0 ;        /* flag whether to use default tmp vols */
  int num_temp_vols;              /* number of tmp drives ava */
  int temp_vol[TMP_VOL_MAX] ;     /* vector of physical tmp drive numbers */
  int get_tempinfo() ;            /* func which reads temp drive datafile */
  int junk ;                      /* temp variable */
/* CDS 9/93 End of new variables */

  /* Look at the environment to see if this is a batch job. */
  if (!(mbatchc = getenv("MBATCH="))) {
    job_id = 0;
    }
  else {
    job_id = atoi(mbatchc);
    sprintf(cleanfile,"Batch_Clean.%05d",job_id);
    /* If the file has not yet been open, then open it and
     * tell it that we want it to delete itself. */
    if (stat(cleanfile,&statjunk)!=0) {
      cleanfilep = fopen(cleanfile,"a");
      if (!cleanfilep) {
        fprintf(stderr,"sequential I/O: couldn't open %s\n",cleanfile);
        }
      else {
        fseek(cleanfilep,(long)0,SEEK_END);
        fprintf(cleanfilep,"%s\n",cleanfile);
        fclose(cleanfilep);
        }
      }
    }

  strcpy(name_format,"%s");

  /* Allocate memory for the unit descriptor. */
  ud = (sequential_t *) malloc(sizeof(sequential_t));
  malloc_check(open_name,ud);

  osi = oldstyleinput();

  /* Find out if we want extra information to be printed about this
   * unit.  name is used as a temporary buffer here. */
  if (osi) {
    strcpy(param,baseparam);
    strcat(param,"verbose");
    if (get_param(param,"%s",name)== -1) ud->verbose = 0;
    else {
      if (!strcmp(name,"on")) ud->verbose = 1;
      else if (!strcmp(name,"yes")) ud->verbose = 1;
      else if (!strcmp(name,"y")) ud->verbose = 1;
      else ud->verbose = 0;
      }
    }
  else {
    if (!read_files_boolean(unit,"VERBOSE",&bool)) ud->verbose = 0;
    else ud->verbose = bool;
    }

  /* See if we want to keep this file when we are done. */
  /* Setting job_id = 0 causes the file to be kept. */
  if (osi) {
    strcpy(param,baseparam);
    strcat(param,"keep");
    if (!(get_param(param,"%s",name)== -1)) {
      if (!strcmp(name,"on")) job_id = 0;
      else if (!strcmp(name,"1")) job_id = 0;
      else if (!strcmp(name,"yes")) job_id = 0;
      else if (!strcmp(name,"y")) job_id = 0;
      }
    }
  else {
    if (read_files_boolean(unit,"KEEP",&bool)) {
      if (bool) job_id = 0;
      }
    }

  /* Set up the block size for this file system. */
  if (osi) {
    strcpy(param,baseparam);
    strcat(param,"blocksize");
    if (get_param(param,"%d",&ud->blocksize) == -1) ud->blocksize = 0;
    }
  else {
    if (!read_files_integer(unit,"BLOCKSIZE",&ud->blocksize)) ud->blocksize = 0;
    }
  if (ud->blocksize == 0) ud->blocksize = 8192;

  /* Find out how many volumes to place the unit across. */
  if (osi) {
    strcpy(param,baseparam);
    strcat(param,"n");
    if (get_param(param,"%d",&ud->n) == -1) {
       ud->n = 1;
       use_default_vols = 1 ;
       }
    }
  else {
    if (!read_files_integer(unit,"NVOLUME",&ud->n)) {
       ud->n = 1;
       use_default_vols = 1 ;
       }
    }
    if (ud->n == 0) { ud->n = 1; use_default_vols = 1; }


  /* Find out how the files are to be named. */
  if (osi) {
    if (get_param("FILES:name",name_format,name) == -1) strcpy(name,"noname");
    }
  else {
    if (!read_files_string(unit,"NAME",name)) strcpy(name,"noname");
    }

/* CDS 9/93 modified code begins here -- check how many temp drives ava */
  got_first_volpath = 1 ;
  if (osi) {
    sprintf(param,"%s%d",baseparam,1);
    if (get_param(param,"%s",volpath) == -1) got_first_volpath = 0 ;
    }
  else {
    sprintf(param,"VOLUME%d",1);
    if (!read_files_string(unit,param,volpath)) got_first_volpath = 0 ;
    }

  /* fprintf(stdout, "sequential_ioopen: Got volpath1 = %s\n", volpath); */

  if (got_first_volpath && use_default_vols) {
     if (strncmp("/tmp",volpath,4)==0) { /*ch defaults for tmp files only*/
        if (!get_tempinfo(&num_temp_vols, temp_vol) ) {
           /* if we couldn't read the datafile, set it this way */
           fprintf(stderr, 
               "sequential_ioopen: can't read default tmp drive file\n") ;
           io_q_abort() ;
           }
        ud->n = num_temp_vols ;
        }
     }
/* CDS 9/93 Finished getting temp drive info and setting num of volumes */

  for (i=0; i<ud->n; i++) {
    if (osi) {
      sprintf(param,"%s%d",baseparam,i);
      if (get_param(param,"%s",volpath) == -1) {
        if (ud->n > 1) no_path_given(open_name);
        else volpath[0] = '\0';
        }
      }
    else {
      sprintf(param,"VOLUME%d",i+1);
      if (!read_files_string(unit,param,volpath)) {
        if (ud->n > 1) no_path_given(open_name);
        else volpath[0] = '\0';
        }
      }

/* CDS 9/93 Now check again if tmp file.  If so, fix vol number           */
   if (use_default_vols && (strncmp("/tmp", volpath, 4)==0)) {
      sscanf(volpath, "/tmp%d%s", &junk, endofvolpath) ;
      sprintf(volpath, "/tmp%d%s", temp_vol[i], endofvolpath) ;
      }
/* CDS 9/93 Ok, that's all the changes except for the new function at end */


    sprintf(path,"%s%s%c%d",volpath,name,'.',unit);
    ud->v[i].path = (char *) malloc(strlen(path)+1);
    malloc_check(open_name,ud->v[i].path);
    strcpy(ud->v[i].path,path);
    /* Any file for with keep is given is kept, otherwise, any
     * file which has a /tmp as the first four characters in the
     * path is keep. */
    if (job_id&&(!strncmp("/tmp",path,4))) {
      cleanfilep = fopen(cleanfile,"a");
      if (!cleanfilep) {
        fprintf(stderr,"sequential I/O: couldn't open %s\n",cleanfile);
        }
      else {
        fseek(cleanfilep,(long)0,SEEK_END);
        fprintf(cleanfilep,"%s\n",path);
        fclose(cleanfilep);
        }
      }
#if defined(DEC)
    if (stat(ud->v[i].path,&statjunk)==0) mode = "r+";
    else mode = "a+";
    ud->v[i].stream = fopen(ud->v[i].path,mode);
    fopen_check(open_name,ud->v[i].path,ud->v[i].stream);
#else
/* #if defined(AIX) || defined(SGI) */
    if (stat(ud->v[i].path,&statjunk)==0) {
      mode = "r+";
      ud->v[i].stream = fopen(ud->v[i].path,mode);
      fopen_check(open_name,ud->v[i].path,ud->v[i].stream);
      }
    else {
      mode = "a+";
      ud->v[i].stream = fopen(ud->v[i].path,mode);
      fopen_check(open_name,ud->v[i].path,ud->v[i].stream);
      fclose(ud->v[i].stream);
      mode = "r+";
      ud->v[i].stream = fopen(ud->v[i].path,mode);
      fopen_check(open_name,ud->v[i].path,ud->v[i].stream);
      }
#endif /* end defined(DEC) */
    }
  ud->next = -1;
  ud->unit = unit;

  if (ud->verbose) {
    fprintf(stderr,"SEQ_IO: opened unit %d {\n",unit);
    fprintf(stderr,"  blocksize = %d\n",ud->blocksize);
    for (i=0; i<ud->n; i++) {
      fprintf(stderr,"  v[%d].path = \"%s\"\n",i,ud->v[i].path);
      }
    fprintf(stderr,"  }\n");
    }
  return(ud);
  }

void
sequential_ioclos(ud,status)
sequential_t *ud;
int status;
{
  int i;

  for (i=0; i<ud->n; i++) {
    fclose(ud->v[i].stream);
    if (status == 4) unlink(ud->v[i].path);
    free(ud->v[i].path);
    ud->v[i].path = NULL;
    }
  if (ud->verbose) {
    fprintf(stderr,"SEQ_IO: closed unit %d {\n",ud->unit);
    fprintf(stderr,"  incount = %u\n",ud->incount);
    fprintf(stderr,"  outcount = %u\n",ud->outcount);
    fprintf(stderr,"  }\n");
    }
  }

void
sequential_iordr(ud,buffer,first,length)
sequential_t *ud;
char *buffer;
unsigned long int length;
PSI_FPTR first;
{
  ud->incount += length;
  sequential_iordwrr("sequential_iordr",IOOP_READ,ud,buffer,first,length);
  }

void
sequential_iowrr(ud,buffer,first,length)
sequential_t *ud;
char *buffer;
unsigned long int length;
PSI_FPTR first;
{
  ud->outcount += length;
  sequential_iordwrr("sequential_iowrr",IOOP_WRITE,ud,buffer,first,length);
  }

void
sequential_iordwrr(caller,ioop,ud,buffer,first,length)
char *caller;
int ioop;
sequential_t *ud;
char *buffer;
unsigned long int length;
PSI_FPTR first;
{
  PSI_FPTR i;
  PSI_FPTR firstbyte, lastbyte;
  PSI_FPTR firstblock, lastblock;
  PSI_FPTR offset;
  PSI_FPTR ncycles;
  PSI_FPTR remainingbytes;
  PSI_FPTR remainingseek;
  PSI_FPTR fullvol;
  PSI_FPTR offset1, offset2;
  PSI_FPTR ibuf;
  PSI_FPTR len;

  firstbyte = first;
  lastbyte = first + length - 1;

  ncycles = firstbyte/(ud->n*ud->blocksize);
  remainingseek = firstbyte -
	ncycles*((PSI_FPTR) ud->n)*((PSI_FPTR) ud->blocksize);
  fullvol = remainingseek/ud->blocksize;

  offset2 = ncycles * ud->blocksize;
  offset1 = offset2 + ud->blocksize;

  /* Seek all volumes to the appropiate positions. */
  if ((ud->next != firstbyte)||(ud->last_ioop != ioop)) {
    for (i=0; i<ud->n; i++) {
      long offset;
      if (i < fullvol) offset = offset1;
      else if (i == fullvol) offset = offset2 + remainingseek%ud->blocksize;
      else offset = offset2;
#if DEBUG
      fprintf(stdout,"seeking volume %d to %ld\n",i,offset);
#endif
      if (fseek(ud->v[i].stream, offset, SEEK_SET)) {
        fprintf(stderr,"%s: fseek: offset = %ld, vol = %d\n",caller,offset,i);
#if defined(AIX)
        fprintf(stderr,"errno = %d\n",errno);
#else
        perror(caller);
#endif
        io_q_abort();
        }
      }
    if (length > 0) ud->last_ioop = ioop;
    }
  ud->next = lastbyte + 1;

  /* Do the io. */
  i = fullvol;
  len = ud->blocksize - remainingseek%ud->blocksize;
  remainingbytes = lastbyte - firstbyte + 1;
#if DEBUG
  fprintf(stdout,"%s: len=%d,remainingbytes=%ld,firstbyte=%ld,lastbyte=%ld,i=%d\n",
          caller,len,remainingbytes,firstbyte,lastbyte,i);
#endif
  ibuf = 0;
  while (remainingbytes > 0) {
    if (len > remainingbytes) len = remainingbytes;
#if DEBUG
    fprintf(stdout,"       len=%d,remainingbytes=%ld,i=%d\n",
          len,remainingbytes,i);
#endif
    if (ioop == IOOP_READ) {
      if (fread(&buffer[ibuf],len,1,ud->v[i].stream)!=1) {
        fprintf(stderr,"%s: unit = %d, len = %d, volume = %d\n",
                caller,ud->unit,len,i);
        if (ferror(ud->v[i].stream)) fread_error(caller);
        }
      }
    else if (ioop == IOOP_WRITE) {
      if (fwrite(&buffer[ibuf],len,1,ud->v[i].stream)!=1) {
        fprintf(stderr,"%s: unit = %d, len = %d, volume = %d\n",
                caller,ud->unit,len,i);
        if (ferror(ud->v[i].stream)) fwrite_error(caller);
        }
      }
    else {
      fprintf(stderr,"%s: illegal ioop = %d\n",caller,ioop);
      io_q_abort();
      }
    i++;
    if (i == ud->n) i=0;
    remainingbytes -= len;
    ibuf += len;
    len = ud->blocksize;
    }
  }

/*
** CDS 9/93
** GET_TEMPINFO : David Sherrill, April 1993
**
** This function will allow PSI to figure out how many temp drives
** to use for sequential io.  The data will be contained in a
** host table file listing each host, the number of temp drives
** it has, and the number labels for each of these drives
**
** Arguments:
**    num_temp_vols = ptr to number of temp vols to use
**    temp_vol      = array of physical temp disks to use
**                      (i.e. {1 2} to use /tmp1 and /tmp2)
**
** Returns: 1 for success, 0 otherwise
*/


#define HOSTNAME_MAX 26

int
get_tempinfo(num_temp_vols, temp_vol)
      int *num_temp_vols, *temp_vol ;
{
FILE *fpi ;                        /* for reading in the host table data */
char hostname[HOSTNAME_MAX] ;      /* name of machine we're running on */
char *hostfile;                    /* filename containing tmp disk info */
char line[MAX_STRING] ;            /* hold line from hostfile */
int found = 0 ;                    /* is host found in data file ? */
char *sptr ;                       /* keep place in input string */
int i, data_in ;

   hostfile = SITEDIR "/tmpdisks.dat" ; 

/* open datafile on hosts' temp disks */
   if ((fpi = fopen(hostfile, "r")) == NULL) {
      fprintf(stderr, "get_tempinfo: couldn't open %s\n", hostfile);
      return(0) ;
      }

   /* get hostname */
   if (gethostname(hostname,HOSTNAME_MAX) == -1) {
      fprintf(stderr, "get_tempinfo: trouble getting hostname\n") ;
      return(0) ;
      }

   /* fprintf(stdout, "get_tempinfo: got hostname = %s\n", hostname) ; */

   /* scan for that hostname in the datafile, get the info */
   while (io_getline(fpi, line) != -1) {
      if (strstr(line, hostname)) { found = 1 ; break ; }
      }
   if (!found) {
      fprintf(stdout, "get_tempinfo: no host %s in datafile\n", hostname) ;
      return(0) ;
      }

   else { /* get the info */
      if ( (sptr = strchr(line, '=')) == NULL) {
         fprintf(stderr, "get_tempinfo: %s has bad format\n", hostfile) ;
         return(0) ;
         }
      sptr++ ;
      while ( (*sptr == ' ') && (*sptr != '\0') ) sptr++ ;

      if (sscanf(sptr, "%d", num_temp_vols) != 1) {
         fprintf(stderr, "get_tempinfo: %s has bad format\n", hostfile) ;
         return(0) ;
         }

      if (*num_temp_vols > TMP_VOL_MAX) {
         fprintf(stderr, "get_tempinfo: %d exceeds %d maximum vols\n",
            *num_temp_vols, TMP_VOL_MAX) ;
         *num_temp_vols = TMP_VOL_MAX ;
         }

      for (i=0; i<(*num_temp_vols); i++) {
         while ( (*sptr != ' ') && (*sptr != '\0') ) sptr++ ;
         while ( (*sptr == ' ') && (*sptr != '\0') ) sptr++ ;
         if (*sptr == '\0') {
            fprintf(stderr, "get_tempinfo: %s bad format\n", hostfile) ;
            return(0) ;
            }
         if (sscanf(sptr, "%d", &data_in) != 1) {
            fprintf(stderr, "get_tempinfo: %s has bad format\n", hostfile) ;
            return(0) ;
            }
         else temp_vol[i] = data_in ;
         }

/*
   fprintf(stdout, "get_tempinfo: got %d vols : ", *num_temp_vols) ;
   for (i=0; i<(*num_temp_vols); i++) {
      fprintf(stdout, "%d ", temp_vol[i]) ;
      }
   fprintf(stdout, "\n") ;
*/

      return(1) ;
      }
}


#if defined(DEC)
char *strstr (s1, s2)    /* for the DEC compiler which lacks strstr() */
      char *s1, *s2 ;
{

   if (*s2 == '\0')
      return ((char *)s1) ;
   for (; (s1 = strchr(s1, *s2)) != NULL; ++s1)
      {
      char *sc1, *sc2 ;

      for (sc1 = s1, sc2 = s2; ; )
         if (*++sc2 == '\0')
            return ((char *)s1) ;
         else if (*++sc1 != *sc2)
            break ;
      }

   return(NULL) ;
}
#endif

/* CDS 9/93 : End modifications */


