/********************************************************************
**
** File:    sat_data.h
**
** Written: Thayaparan Thanabalasingham
**
** Last Modified:    30. 06. 1998 - 01:36:00
**
********************************************************************/


#define LINE 80
#define NAME 12


/* Structure that holds three lines of nasa format data */
typedef struct Line_struct {

   char line0[LINE];
   char line1[LINE];
   char line2[LINE];

} LineStruct;


/* Structure that holds satellite data */
typedef struct SatData_struct {
   
   char satName[NAME];
   
   unsigned long satNumber;
   unsigned long elementSet;
   unsigned long epochRev;
   
   double epochTime;
   double epochYear;
   double epochRaan; 
   double epochMeanMotion;
   double epochMeanAnomaly;
   double epochArgOfPerigee;
   double decayRate;
   double inclination;
   double eccentricity;

} SatDataStruct;



int getSatData (const char *satName, SatDataStruct *sdp);
void printData (const SatDataStruct *sdp);

extern FILE *fsat;

