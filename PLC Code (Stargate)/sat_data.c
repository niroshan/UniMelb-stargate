      /*****************************************************************
      **                                                              **
      ** File: sat_data.c                                             **
      **                                                              **
      ** Written: Thayaparan Thanabalasingham                         **
      **                                                              **
      ** Last Modified: 04. 07. 1998 - 19:45:00                       **
      **                                                              **
      *****************************************************************/


   /*************************************************************************
   **                                                                      **
   **      Extracting satellite data from NASA format input lines          **
   **      ======================================================          **
   **                                                                      **
   ** The format consists of groups of 3 lines: One line containing        **
   ** the satellite name, followed by the standard Two-Line Orbital        **
   ** Element Set Format identical to that used by NASA and NORAD.         **
   **                                                                      **
   ** Line 1                                                               **
   ** ======                                                               **
   **                                                                      **
   ** Line 1 is a eleven-character name.                                   **
   **                                                                      **
   ** -------------------------------------------------------------------- **
   **                                                                      **
   ** Line 2                                                               **
   ** ======                                                               **
   **                                                                      **
   ** Column       Description                                             **
   **                                                                      **
   **  01-01       Line Number of Element Data                             **
   **  03-07       Satellite Number                                        **
   **  10-17       International Designator                                **
   **  19-20       Epoch Year (Last two digits of year)                    **
   **  21-32       Epoch (Julian Day and fractional portion of the day)    **
   **  34-43       First Time Derivative of the Mean Motion divided by 2.  **
   **           or Ballistic Coefficient. (Depending of ephemeris type)    **
   **  45-52       Second Time Derivative of Mean Motion divided by 6.     **
   **              (Blank if N/A)                                          **
   **  54-61       BSTAR drag term if GP4 general perturbation is theory   **
   **              was used. Otherwise, radiation pressure coefficient.    **
   **  63-63       Ephemeris type                                          **
   **  65-68       Element Number                                          **
   **  69-69       Check Sum (Modulo 10)                                   **
   **              (Letters, blanks, periods = 0; minus sign = 1;          **
   **              plus sign = 0)                                          **
   **                                                                      **
   **                                                                      **
   ** All other colums are blank or fixed.                                 **
   **                                                                      **
   ** -------------------------------------------------------------------- **
   **                                                                      **
   ** Line 3                                                               **
   ** ======                                                               **
   **                                                                      **
   ** Column       Description                                             **
   **                                                                      **
   **  01-01       Line Number of Element Data                             **
   **  03-07       Satellite Number                                        **
   **  09-16       Inclination [Degrees]                                   **
   **  18-25       Right Ascension of the Ascending Node [Degrees]         **
   **  27-33       Eccentricity (decimal points assumed)                   **
   **  35-42       Argument of Perigee [Degrees]                           **
   **  44-51       Mean Anomaly [Degrees]                                  **
   **  53-63       Mean Motion [Recs per day]                              **
   **  64-68       Revolution Number at epoch [Revs]                       **
   **  69-69       Check Sum (Modulo 10)                                   **
   **                                                                      **
   **                                                                      **
   ** All other colums are blank or fixed.                                 **
   **                                                                      **
   *************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "sat_data.h"
#include "sat_track.h"
#include "utilities.h"
#include "general.h"


           /* Functions local to this module */

/* ########## - Just to shut up GNUC ############# */
int strncasecmp (const char *p1, const char *p2, size_t size); 

static int findData (const char *sat, LineStruct *lp);
static int checkData (const LineStruct *lp);
static int testCheckSum (const char *line); 
static double getValue (const char line[], int start, int end);


FILE *fsat;



/*****************************************************************
**
** getSatData ().
**
**      Extracts the satellite data for a given satellite.
**
*****************************************************************/

int getSatData (const char *satName, SatDataStruct *sdp)
{
   LineStruct l;


   /* Error cases: Satellite data not found in the file OR Check sum error */
   if ((findData (satName, &l) == 0) || (checkData (&l) == 0))
      return 0;

   /* Get the Satellite Parameters */
   strcpy (sdp->satName, l.line0);
   sdp->satNumber  = (long) getValue (l.line1,  2,  6);
   sdp->elementSet = (long) getValue (l.line1, 64, 67);
   sdp->epochRev   = (long) getValue (l.line2, 63, 67);

   sdp->epochYear         = getValue (l.line1, 18, 19);
   sdp->epochTime         = getValue (l.line1, 20, 31);
   sdp->decayRate         = getValue (l.line1, 33, 42);
   sdp->inclination       = getValue (l.line2,  8, 15) * RADIANS_PER_DEGREE;
   sdp->epochRaan         = getValue (l.line2, 17, 24) * RADIANS_PER_DEGREE;
   sdp->eccentricity      = getValue (l.line2, 26, 32);
   sdp->epochArgOfPerigee = getValue (l.line2, 34, 41) * RADIANS_PER_DEGREE;
   sdp->epochMeanAnomaly  = getValue (l.line2, 43, 50) * RADIANS_PER_DEGREE;
   sdp->epochMeanMotion   = getValue (l.line2, 52, 62);

   /* Number of days from Reference Time (1900 Jan 01) */
   sdp->epochTime    += getDayNum (sdp->epochYear, 1, 0);
   /* decimal point - see comments on NASA format */
   sdp->eccentricity *= 1.0e-7;

   return 1;
}


/*************************************************************************
**
** findData ().
**
** Finding the "two-line" format satellite-data from the database.
**
** Assumes that the file does not have any blank lines or comments.
** ("MAY WANT TO CHANGE THAT LATER." 04. 07. 98 19:35:00)
**
*************************************************************************/

int findData (const char *satName, LineStruct *lp)
{
   char line[LINE];

        /*Look for the satellite name in the database.*/
   fgets (line, LINE, fsat);

   while (strncasecmp (line, satName, strlen (satName)) != 0) {
      /* skip next two lines */
      fgets (line, LINE, fsat);
      fgets (line, LINE, fsat);

      /* Read the next satellite name from the file */
      if ((fgets (line, LINE, fsat) == NULL)) {
         /* Couldn't find the satellite in the data base. */
         printf ("\n\tSatellite \"%s\" not found.\n\n", satName);
         rewind (fsat);
         return 0;
      }
   }

   strcpy (lp->line0, satName);
   fgets (lp->line1, LINE, fsat);
   fgets (lp->line2, LINE, fsat);

   rewind (fsat);

   return 1;
}


/*****************************************************************
**
** checkData ().
**
**   Testing the "two line" format satellite-data for errors.
**
**   Checks include,
**      i) Line numbers of element data.
**     ii) Checksum tests.
**
*****************************************************************/

int checkData (const LineStruct *lp)
{
             /* Checking line numbers. */
   if (strncmp (lp->line1, "1 ", 2) != 0) {
      fprintf (stdout, "Line 1 is not present\n");
      return 0;
   }
   if (strncmp (lp->line2, "2 ", 2) != 0) {
      fprintf (stdout, "Line 2 is not present\n");
      return 0;
   }

                /* Check sum tests. */
   if (testCheckSum (lp->line1) == 0) {
      fprintf (stdout, "Check sum error on line 1\n");
      return 0;
   }
   if (testCheckSum (lp->line2) == 0) {
      fprintf (stdout, "Check sum error on line 2\n"); 
      return 0;
   }
   
   return 1;
}
   

/*******************************************************************
**
** testCheckSum ().
**
**           Do a check sum test on the given line.
**
** - Last digit (69th) of the line is a modulo-10 checksum digit.
** - Letters, blanks, periods, plus-sign = 0; minus sign = 1.
**
*******************************************************************/

static int testCheckSum (const char *line) 
{
   int i;
   int sum = 0;
   int checkSum;

   checkSum = (line[68] - '0');

   for (i = 0; i < 68; i++) {
      if (line[i] == '-')
         sum += 1; 
      else if (isdigit (line[i]))
         sum += (line[i] - '0');
   }

   return ((sum % 10) == checkSum);
}


/*****************************************************************
**
** getValue ().
**
**         Converts the given string to a number (double).
**
*****************************************************************/

static double getValue (const char line[], int start, int end)
{
   int len;
   char string[LINE];
   double value;

   len = end - start + 1;
   
   strncpy (string, line + start, len);
   string[len] = '\0';

   value = atof (string);

   return value;
}


/*****************************************************************
**
** printData ().
**
**  Printing the Satellite Data in AMSAT Format (For Testing).
**
*****************************************************************/

void printData (const SatDataStruct *sdp) 
{
   fprintf (stdout, "\n\n");

   fprintf (stdout,
         "Satellite:\t %s\n"       , sdp->satName);
   fprintf (stdout,
         "Catalog number:\t %lu\n" , sdp->satNumber);
   fprintf (stdout,
         "Element set:\t %lu\n"    , sdp->elementSet);
   fprintf (stdout, 
         "Epoch time:\t %.8f\n"    , sdp->epochTime);
   fprintf (stdout, 
         "Inclination:\t %.4f\n"   , sdp->inclination * DEGREES_PER_RADIAN);
   fprintf (stdout, 
         "RA of node:\t %.4f\n"    , sdp->epochRaan * DEGREES_PER_RADIAN);
   fprintf (stdout, 
         "Eccentricity:\t %.7f\n"  , sdp->eccentricity);
   fprintf (stdout, 
         "Arg of perigee:\t %.4f\n"
         , sdp->epochArgOfPerigee * DEGREES_PER_RADIAN);
   fprintf (stdout, 
         "Mean anomaly:\t %.4f\n"  
         , sdp->epochMeanAnomaly * DEGREES_PER_RADIAN);
   fprintf (stdout, 
         "Mean motion:\t %.8f\n"   , sdp->epochMeanMotion);
   fprintf (stdout, 
         "Decay rate:\t %.8f\n"    , sdp->decayRate);
   fprintf (stdout, 
         "Epoch rev:\t %lu\n"      , sdp->epochRev);

   fprintf (stdout, "\n\n");
}

