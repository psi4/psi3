/*!
** \file get_numvols.c
*/
 
#include <stdio.h>
#include <string.h>
#include <libipv1/ip_lib.h>
#include "psio.h"

int psio_get_tempinfo(ULI *num_temp_vols);

/*!
** PSIO_GET_NUMVOLS(): Get the number of volumes that file number 'unit'
** is split across.
*/
ULI psio_get_numvols(ULI unit)
{
  ULI num;
  int errcod;
  char ip_token[PSIO_MAXSTR];
  char *gprgid();

  sprintf(ip_token,":%s:FILES:FILE%u:NVOLUME",gprgid(),unit);
  errcod = ip_data(ip_token,"%u",&num,0);
  if(errcod == IPE_OK) return(num);

  sprintf(ip_token,":%s:FILES:DEFAULT:NVOLUME",gprgid());
  errcod = ip_data(ip_token,"%u",&num,0);
  if(errcod == IPE_OK) return(num);

  sprintf(ip_token,":DEFAULT:FILES:FILE%u:NVOLUME",unit);
  errcod = ip_data(ip_token,"%u",&num,0);
  if(errcod == IPE_OK) return(num);

  sprintf(ip_token,":DEFAULT:FILES:DEFAULT:NVOLUME");
  errcod = ip_data(ip_token,"%u",&num,0);
  if(errcod == IPE_OK) return(num);

  if(!psio_get_tempinfo(&num)) return(0);

  return(num);
}


/*!
** PSIO_GET_NUMVOLS_DEFAULT(): Get the number of volumes that file 
** number 'unit' is split across.
*/
ULI psio_get_numvols_default(void)
{
  ULI num;
  int errcod;
  char ip_token[PSIO_MAXSTR];

  sprintf(ip_token,":DEFAULT:FILES:DEFAULT:NVOLUME");
  errcod = ip_data(ip_token,"%u",&num,0);
  if(errcod == IPE_OK) return(num);

  if(!psio_get_tempinfo(&num)) return(0);

  return(num);
}



/*!
** PSIO_GET_TEMPINFO : David Sherrill, April 1993
**
** This function will allow PSI to figure out how many temp drives
** to use for sequential io.  The data will be contained in a
** host table file listing each host, the number of temp drives
** it has, and the number labels for each of these drives
**
** \param num_temp_vols = ptr to number of temp vols to use
**  
** Returns: 1 for success, 0 otherwise
*/

#define HOSTNAME_MAX 26

int psio_get_tempinfo(ULI *num_temp_vols)
{
   FILE *fpi ;                        /* for reading in the host table data */
   char hostname[HOSTNAME_MAX] ;      /* name of machine we're running on */
   char *hostfile;                    /* filename containing tmp disk info */
   char line[PSIO_MAXSTR] ;         /* hold line from hostfile */
   int found = 0 ;                    /* is host found in data file ? */
   char *sptr ;                       /* keep place in input string */
   int i, data_in ;

   hostfile = SITEDIR "/tmpdisks.dat" ;

   /* open data file on hosts' temp disks */
   fpi = fopen(hostfile, "r");

   /* open datafile on hosts' temp disks */
   if (fpi == NULL) {
      fprintf(stderr, "get_tempinfo: couldn't open %s\n", hostfile);
      fclose(fpi);
      return(0) ;
      }

   /* get hostname */
   if (gethostname(hostname,HOSTNAME_MAX) == -1) {
      fprintf(stderr, "get_tempinfo: trouble getting hostname\n") ;
      fclose(fpi);
      return(0) ;
      }

   /* fprintf(stdout, "get_tempinfo: got hostname = %s\n", hostname) ; */

   /* scan for that hostname in the datafile, get the info */
   while (io_getline(fpi, line) != -1) {
      if (strstr(line, hostname)) { found = 1 ; break ; }
      }
   if (!found) {
      fprintf(stdout, "get_tempinfo: no host %s in datafile\n", hostname) ;
      fclose(fpi);
      return(0) ;
      }

   else { /* get the info */
      if ( (sptr = strchr(line, '=')) == NULL) {
         fprintf(stderr, "get_tempinfo: %s has bad format\n", hostfile) ;
	 fclose(fpi);
         return(0) ;
         }
      sptr++ ;
      while ( (*sptr == ' ') && (*sptr != '\0') ) sptr++ ;

      if (sscanf(sptr, "%u", num_temp_vols) != 1) {
         fprintf(stderr, "get_tempinfo: %s has bad format\n", hostfile) ;
	 fclose(fpi);
         return(0) ;
         }

      if (*num_temp_vols > PSIO_MAXVOL) {
         fprintf(stderr, "get_tempinfo: %u exceeds %u maximum vols\n",
            *num_temp_vols, PSIO_MAXVOL) ;
         *num_temp_vols = PSIO_MAXVOL ;
         }

      fclose(fpi);
      return(1) ;
      }
}
