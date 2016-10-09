
/**************************************************************
**
** File: "utilities.c"
**
** Written: Thayaparan Thanabalasingham
**
** Last Modified: 12.07.98 17:56:00
**
**
**************************************************************/


#include <stdio.h>
#include <math.h>
#include <time.h>
#include "utilities.h"
#include "general.h"


/* Day numbers for start of each month from Jan 01 */
int monthDays [] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};


/*****************************************************************
**
** getDayNum ().
**
** Finding the Day Number for a given date. January 1 of the
** reference year is day 0.
**
** January 1, 1900 is day 0
** valid from 1901 through 2099.
** (Between 1901 and 2099, year divisible by 4 is a leap year).
**
*****************************************************************/

long getDayNum (int year, int month, int day)
{
   long result;

   /* Heuristic to allow 4 or 2 digit year specifications */
   if (year < 50)
      year += 2000;
   else if (year < 100)
      year += 1900;

   result = ((((long) year - 1901) * 1461) >> 2)
          + monthDays[month - 1] + day + 365;

   if (year % 4 == 0 && month > 2)
      result++;

   return result;
}

   
/*************************************************************************
**
** getTime ().
**
** Get the current time (GMT).
**
*************************************************************************/
/*
void getTime (int *day, int *month, int *year, int *hour, int *min, int *sec)
{
    time_t timeofday;
    struct tm *tm;

    time (&timeofday);                
    tm = gmtime (&timeofday);        

    *day   = tm->tm_mday + 1;
    *month = tm->tm_mon + 1;
    *year  = tm->tm_year;
    *hour  = tm->tm_hour;
    *min   = tm->tm_min;
    *sec   = tm->tm_sec;

    return;
}
*/

double getTime ()
{
   int year, month, day, hour, min, sec; 
   time_t timeofday;
   struct tm *tm;
   double value;

   time (&timeofday);                  /* get the UNIX time              */
   tm = gmtime (&timeofday);           /* Greenwich Mean Time            */ 

   day   = tm->tm_mday + 1;
   month = tm->tm_mon + 1;
   year  = tm->tm_year;
   hour  = tm->tm_hour;
   min   = tm->tm_min;
   sec   = tm->tm_sec;

   value = getDayNum (year, month, day)
      + (((sec / 60.0 + min) / 60.0 + hour) / 24.0);

   return value;
}

/**************************************************************************
**
** getDate ().
**
** Finding the date from the day number (number of days from 1900 Jan 01).
**
**************************************************************************/

void getDate (long dayNum, int *year, int *month, int *day)
{
   int m, l;
   long y;

   y  = 4 * dayNum;
   y /= 1461;

   dayNum =  dayNum - 365 - (((y - 1) * 1461) >> 2);

   l = 0;

   if (y % 4 == 0 && dayNum > monthDays[2])
      l = 1;

   m = 1;

   while (dayNum > monthDays[m] + l)
      m++;

   dayNum -= (monthDays[m - 1]);

   if (m > 2)
      dayNum -= l;

   *year  = y + 1900;
   *month = m;
   *day   = dayNum;
}  


/*********************************************************************
**
** printTime ().
**
** Print the date as a string.
**
**********************************************************************/

void printTime (double t)
{
   long d;
   int year, month, day, hour, min, sec; 
   double frac;

   d    = floor (t);
   frac = t - d;
  
   frac *= 24.0;
   hour  = floor (frac);
   frac  = (frac - hour) * 60;
   min   = floor (frac);
   sec   = floor ((frac - min) * 60);
   
   getDate (d, &year, &month, &day); 

   fprintf (stdout, "%d.%02d.%02d %02d:%02d:%02d ",
         year, month, day, hour, min, sec);
}


/* Dummy fn */
void enable ()
{

}


/* Dummy fn */
void disable ()
{

}


void debug_going_down ()
{
   debug_spaces += 2;
}

void debug_going_up ()
{
   debug_spaces -= 2;
}

void print_spaces ()
{
   int i;
   for (i=0; i < debug_spaces ; i++)
      printf (" ");
     /* putchar  (' '); */
}
