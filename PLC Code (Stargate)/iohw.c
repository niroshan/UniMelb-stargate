
/**************************************************************************
**
** file: iohw.c
**
***************************************************************************/


#include <stdlib.h>
#include "iohw.h"


/****************************************************************************
**
** OutputVolts ().
**
****************************************************************************/

int OutputVolts (double voltage) 
{
   double min = -5;
   double max = 5;

   int volt;  /* convert float point voltage */

   voltage=(voltage - min) / (max - min) * 4096;

   if (voltage < 0.0)
      voltage = 0.0;    /* 0 <= volt <= 4095 */
   else
      if (voltage > 4095.0)
         voltage = 4095.0;

   volt = (int) voltage;

   return volt;

} /* end of OutputVolts */


/*************************************************************************
**
** convert_encoder ().
**
************************************************************************/

double convert_encoder (unsigned char type, int value)
{
  double x = 0.0;

  if (type == AZ) {
     x = (double) value;
  }

  if (type == EL) {
     x = (double) value;
  }

  return x;
}

