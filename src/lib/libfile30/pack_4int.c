/* pack_4int()   Packs a matrix of small integers (normally, symmetry information)
**               into a sequential array of integers.
**
**  FORMAT: source is converted into a packed form. Each element of a row is
**  represented by a byte, therefore if nirreps <= 4 - each row is packed   
**  into an integer word, else if nirreps == 8 - into two integer words.
**   Here's how it works :
**   | source[3] | source[2] | source[1] | source[0] |                 
**   leftmost byte                       rightmost byte
**			
**
**  arguments: int **source - integer matrix num_rows by nirreps
**             int  *dest   - the destination array of
**                            num_rows*((nirreps == 8) ? 2 : 1) integers
**             int num_rows - number of rows in source
**             int nirreps  - number of columns in source
**
**  returns: nothing	
*/

void pack_4int(int **source, int *dest, int num_rows, int nirreps)
{
  int i;

  switch (nirreps) {
    case 8: bzero((char *) (dest+num_rows),sizeof(int)*num_rows);
    case 4:
    case 2:
    case 1: bzero((char *) dest,sizeof(int)*num_rows);
  }
  for(i=0;i<num_rows;i++)
    switch (nirreps) {
      case 8: dest[i+num_rows] += source[i][4];
              dest[i+num_rows] += (source[i][5] << 8);
              dest[i+num_rows] += (source[i][6] << 16);
              dest[i+num_rows] += (source[i][7] << 24);
      case 4: dest[i] += (source[i][2] << 16);
              dest[i] += (source[i][3] << 24);
      case 2: dest[i] += (source[i][1] << 8);
      case 1: dest[i] += source[i][0];
    }

}
