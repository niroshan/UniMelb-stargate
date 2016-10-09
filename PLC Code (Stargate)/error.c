#include <stdio.h>
#include "error.h"

/***************************************************************************
**
**                      Thaya 29.07.98
**
***************************************************************************/

void error_set(const char* msg , int code) 
{
   printf ("%s : Error no = %d\n", msg, code);
}
