/* converts integer file pointer to sector pointer */

#include "iomrparam.h"

int i2sec(PSI_FPTR n)
    {
      int num;

      num = (int) ((n+4096)/4096);
      return(num);
    }
