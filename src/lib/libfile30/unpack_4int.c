/* unpack_4int() Unpacks a sequential array of integers into a
**               matrix of small integers (normally, symmetry information)
**
**  FORMAT:  source is converted from its packed form. Each element of a row of dest is
**  represented by a byte, therefore if nirreps <= 4 - each row is packed   
**  into an integer word, else if nirreps == 8 - into two integer words.
**   Here's how it works :
**   | dest[3] | dest[2] | dest[1] | dest[0] |                 
**   leftmost byte                rightmost byte
**			
**
**  arguments: int **dest  - integer matrix num_rows by nirreps
**             int *source - the source array of
**                            num_rows*((nirreps == 8) ? 2 : 1) integers
**             int num_rows - number of rows in dest
**             int nirreps  - number of columns in dest
**
**  returns: nothing	
*/

void unpack_4int(int *source, int **dest, int num_rows, int nirr)
{
  int i;

  for(i=0;i<num_rows;i++)
    switch (nirr) {
      case 8: dest[i][4] = source[i+num_rows] & 255;                   
              dest[i][5] = (source[i+num_rows] >> 8) & 255;
              dest[i][6] = (source[i+num_rows] >> 16) & 255;
              dest[i][7] = (source[i+num_rows] >> 24) & 255;
      case 4: dest[i][2] = (source[i] >> 16) & 255;
              dest[i][3] = (source[i] >> 24) & 255;
      case 2: dest[i][1] = (source[i] >> 8) & 255;
      case 1: dest[i][0] = source[i] & 255;
    }

}
