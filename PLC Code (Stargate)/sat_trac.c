/*********************************************************
**
** File: "sat_track.c"
**
** Written: Thayaparan Thanabalasingham
**
** Last Modified: 29. 07. 98 - 20:54:00
**
*********************************************************/


/*
** Algorithm is obtained from 'n3emo' program.
*/


#include <stdio.h>
#include <math.h>
#include <string.h>

#include "sat_data.h"
#include "sat_track.h"
#include "utilities.h"
#include "general.h"


#define TWOPI            (PI*2)
#define TWOTHIRD         (2.0/3)
#define MINUTESPERDAY    (24 * 60.0)
#define EPSILON          (RADIANS_PER_DEGREE / 3600)     /* 1 arc second */

#define EARTHRADIUS      6378.16              /* Kilometers, at equator */
#define EARTHFLAT        (1 / 298.25)         /* Earth Flattening Coeff. */
#define SIDEREALSOLAR    1.0027379093

#define square(x)        ((x) * (x))


/* Function prototypes for functions local to this module */

/* ########## - Just to shut up GNUC */
int strncasecmp (const char *p1, const char *p2, size_t size); 

static void getPrecession (double a, double e, double i);
static double solveTimeKepler (register double meanAnomaly, 
      register double e, register double a,
      register double *radius);
static void getSatPosition (double inclination, double epochRaan, 
      double epochArgPerigee, double epochTime, double semiMajorAxis,
      double time, double trueAnomaly, double radius,
      double *satX, double *satY, double *satZ);
static void getSitePosition (double latitude, double longitude,
      double altitude, double time,
      double *siteX, double *siteY, double *siteZ,
      mat3x3 siteMatrix);
static void getTopocentric (double satX, double satY, double satZ,
      double siteX, double siteY, double siteZ,
      mat3x3 siteMatrix,
      double *x, double *y, double *z);
static void getAzEl (double x, double y, double z,
      double *azimuth, double *elevation);
static double calcLocalSideReal (double time, double longitude);


/* Used to avoid unneccesary recomputation */
char prevSatName [12] = "";
double oldSiteLat  = -100000;  
double oldSiteLong = -100000;

double raanPrecession;
double perigeePrecession;
double referenceOrbit;

SatDataStruct satData;


/**************************************************************************
**
** satTrack ().
**
** Calculates azimath and elevation of a satellite at a given time.
**
** Calculation involves following steps:
** - calculating orbital plane parameters true anomaly and radius
** - calculating satellite position in right ascension declination system
** - calculating site position in right ascension declination system
** - tranformation from right ascention system to topocentric (observer
**   centered) system
** - calculating azimath and elevation from topocentric coordinates
**
**************************************************************************/

void satTrack (const char satName[], double currentTime,
      const SiteStruct *site, double *azimuth, double *elevation) 
{
   double currentMotion;
   double period;
   double semiMajorAxis;
   
   double averageMotion;
   double currentOrbit;
   double meanAnomaly;
   
   double trueAnomaly;
   double radius;
   
   double satX, satY, satZ;
   double siteX, siteY, siteZ;
   double x, y, z;

   mat3x3 siteMatrix;

   /* This calculation is needed only once for one tracking sequence */
   if (strncasecmp (satName, prevSatName, 11) != 0) {
      getSatData (satName, &satData);

      /* Kepler's 3rd Law: Calculate semimajor axis at epoch time */
      period = MINUTESPERDAY / satData.epochMeanMotion; /* in minutes */
      semiMajorAxis = 331.25 * pow (period, TWOTHIRD);

      /* Compensate for precession of equox vernox 00000 */
      getPrecession (semiMajorAxis, satData.eccentricity, satData.inclination);

      /* Orbit at epoch time */
      referenceOrbit = satData.epochMeanAnomaly / TWOPI + satData.epochRev;
      
      strcpy (prevSatName, satName); 
   }

   /*
   ** Calculate semimajor axis
   ** Kepler's 3rd Law: semimajor axis (a) vs period (T).
   ** a^3 = K * T^2.
   ** a in kilometers
   ** T in minutes and
   ** K = 331.25^3 km^3 / min^2.
   */

   currentMotion = satData.epochMeanMotion 
                 + (currentTime - satData.epochTime) * satData.decayRate * 2;
   period        = MINUTESPERDAY / currentMotion; /* in minutes */
   semiMajorAxis = 331.25 * pow (period, TWOTHIRD);
                 
   /* Calculate Mean Anomaly */
   averageMotion = satData.epochMeanMotion 
                 + (currentTime - satData.epochTime) * satData.decayRate;
   currentOrbit = referenceOrbit 
                + (currentTime - satData.epochTime) * averageMotion;
   meanAnomaly  = (currentOrbit - floor (currentOrbit)) * TWOPI;

   /* Get orbital plane parameters: true anomaly and radius */
   trueAnomaly = solveTimeKepler (meanAnomaly, satData.eccentricity,
         semiMajorAxis, &radius);

   /* Get the satellite position in right ascension declination system */
   getSatPosition (satData.inclination, satData.epochRaan,
         satData.epochArgOfPerigee, satData.epochTime, semiMajorAxis,
         currentTime, trueAnomaly, radius, &satX, &satY, &satZ);

   /* Get the site position in right ascension declination system */
   getSitePosition  (site->latitude, site->longitude, site->altitude, 
         currentTime, &siteX, &siteY, &siteZ, siteMatrix);

   /* Get the satellite position in the topocentric coordinate system */
   getTopocentric (satX, satY, satZ, siteX, siteY, siteZ,
         siteMatrix, &x, &y, &z);
   
   /* Calculate azimath and elevation */
   getAzEl (x, y, z, azimuth, elevation);
}

 
/*************************************************************************
**
** getPrececession ().
**
** Calculated once for each satellite.
**
** a - semimajor axis
** i - inclination
** e - eccentricity
**
************************************************************************/

void getPrecession (double a, double e, double i)
{
   raanPrecession    = 9.95 * pow (EARTHRADIUS / a, 3.5) * cos (i)
                     / square (1 - square (e)) * RADIANS_PER_DEGREE;

   perigeePrecession = 4.97 * pow (EARTHRADIUS / a, 3.5)
                        * (5 * square (cos (i)) - 1)
                        / square (1 - square (e))
                        * RADIANS_PER_DEGREE;
}


 
/*************************************************************************
**
** solveTimeKepler ().
**
** Calculates orbital plane parameters - true anomaly and radius.
**
** Procedure:
** ---------
** Solve Time-Kepler's equation: M = E - e * sin (E).
**      where,
**      e - eccentricity of orbit's ellipse
**      M - mean anomaly
**      E - eccentric anomaly
**
** M and e are known. Use Newton's procedure with first approximation to
** find eccentric anomaly E.
**
** Once, E is known, and semimajor axis is given orbital plane parameters
** can be calculated.
**
** Inputs:
**    meanAnomaly, eccentricity e and semimajoraxis a.
**
** Outputs:
**    trueAnomaly and radius.
**
*************************************************************************/

double solveTimeKepler (register double meanAnomaly, register double e,
      register double a, register double *radius)
{
   register double E;              /* Eccentric Anomaly */
   register double error;
   register double trueAnomaly;

   E = meanAnomaly ;  /* Initial guess */

   do {
      error = (E - e * sin (E) - meanAnomaly) / (1 - e * cos (E));
      E -= error;
   } while (fabs (error) >= EPSILON);

   if (fabs (E - PI) < EPSILON)
      trueAnomaly = PI;
   else
      trueAnomaly = 2 * atan (sqrt ((1 + e) / (1 - e)) * tan (E / 2));

   if (trueAnomaly < 0)
      trueAnomaly += TWOPI;

   (*radius) = a * (1 - square (e)) / (1 + e * cos (trueAnomaly));

   return trueAnomaly;
}
 

/************************************************************************
**
** getSatPosition ().
**
** Compute the satellite postion in the RA based coordinate system 
**
************************************************************************/

void getSatPosition (double inclination, double epochRaan, 
      double epochArgPerigee, double epochTime, double semiMajorAxis,
      double time, double trueAnomaly, double radius,
      double *satX, double *satY, double *satZ)
{
   double raan, argPerigee;
   double xw, yw;	/* In orbital plane */
   double px, qx, py, qy, pz, qz;	/* Escobal transformation 31 */
   double cosAP, sinAP;
   double cosRaan, sinRaan, cosInc, sinInc;

   xw = radius * cos (trueAnomaly);
   yw = radius * sin (trueAnomaly);

   argPerigee = epochArgPerigee + (time - epochTime) * perigeePrecession;
   
   raan       = epochRaan - (time - epochTime) * raanPrecession;

   cosRaan = cos (raan);
   sinRaan = sin (raan);
   cosAP = cos (argPerigee);
   sinAP = sin (argPerigee);
   cosInc = cos (inclination);
   sinInc = sin (inclination);

   px =  cosAP * cosRaan - sinAP * sinRaan * cosInc;
   py =  cosAP * sinRaan + sinAP * cosRaan * cosInc;
   pz =  sinAP * sinInc;
   qx = -sinAP * cosRaan - cosAP * sinRaan * cosInc;
   qy = -sinAP * sinRaan + cosAP * cosRaan * cosInc;
   qz =  cosAP * sinInc;

   (*satX) = px * xw + qx * yw;	/* Escobal, transformation #31 */
   (*satY) = py * xw + qy * yw;
   (*satZ) = pz * xw + qz * yw;
}


/*************************************************************************
**
** getSitePosition ()
**
** Compute the site position in the RA based coordinate system. 
** SiteMatrix is set to a matrix which is used by getTopoCentric ()
** to convert geocentric coordinates to topocentric (observer-centered)
** coordinates.
**
*************************************************************************/

void getSitePosition (double latitude, double longitude, double altitude,
      double time,
      double *siteX, double *siteY, double *siteZ, mat3x3 siteMatrix)
{
   double geodeticLat;   /* Geodetic latitude */
   static double G1 ,G2; /* Used to correct for flattening of the Earth */
   static double cosGL, sinGL;

   double localSideReal; /* local side real time */
   double cosLSR, sinLSR;

   if ((latitude != oldSiteLat) || (longitude != oldSiteLong)) {
      oldSiteLat  = latitude;
      oldSiteLong = longitude;
      geodeticLat = atan (1 / (1 - square (EARTHFLAT)) * tan (latitude));

      cosGL = cos (geodeticLat);
      sinGL = sin (geodeticLat);

      G1 = EARTHRADIUS / 
         (sqrt (1 - (2 * EARTHFLAT - square (EARTHFLAT)) * square (sinGL)));
      G2 = G1 * square (1 - EARTHFLAT);
      G1 += altitude;
      G2 += altitude;
   }

   localSideReal = calcLocalSideReal (time, longitude);

   cosLSR = cos (localSideReal);
   sinLSR = sin (localSideReal);

   (*siteX) = G1 * cosGL * cosLSR;
   (*siteY) = G1 * cosGL * sinLSR;
   (*siteZ) = G2 * sinGL;

   siteMatrix[0][0] = sinGL * cosLSR;
   siteMatrix[0][1] = sinGL * sinLSR;
   siteMatrix[0][2] = -cosGL;
   siteMatrix[1][0] = -sinLSR;
   siteMatrix[1][1] = cosLSR;
   siteMatrix[1][2] = 0.0;
   siteMatrix[2][0] = cosLSR * cosGL;
   siteMatrix[2][1] = sinLSR * cosGL;
   siteMatrix[2][2] = sinGL;
}


/*************************************************************************
** 
** getTopoCentric ().
**
** Convert from geocentric RA based coordinates to topocentric
** (observer centered) coordinates.
**
*************************************************************************/

void getTopocentric (double satX, double satY, double satZ,
      double siteX, double siteY, double siteZ,
      mat3x3 siteMatrix, double *x, double *y, double *z)
{
   satX -= siteX;
   satY -= siteY;
   satZ -= siteZ;

   (*x) = siteMatrix[0][0] * satX 
        + siteMatrix[0][1] * satY
        + siteMatrix[0][2] * satZ; 
   
   (*y) = siteMatrix[1][0] * satX 
        + siteMatrix[1][1] * satY
        + siteMatrix[1][2] * satZ; 
   
   (*z) = siteMatrix[2][0] * satX
        + siteMatrix[2][1] * satY
        + siteMatrix[2][2] * satZ; 
}


/*************************************************************************
**
** getAzEl ().
**
** Given topocentric coordinates, calculates the azimath and elevation.
**
*************************************************************************/

void getAzEl (double x, double y, double z, double *azimuth, double *elevation)
{

   (*elevation) = atan (z / sqrt (square (x) + square (y)));

   (*azimuth)   = PI - atan2 (y, x);

   if ((*azimuth) < 0)
      (*azimuth) += PI;
}


/************************************************************************
** 
** calcLocalSideReal ().
**
** Calculating the sidereal time.
**
************************************************************************/

double calcLocalSideReal (double time, double longitude)
{
   double Tu;
   double sideRealG0;   /* Greenwich sidereal time at 0 hrs UT */
   double sideRealG;    /* Greenwich sidereal time at the given time */
   double localSideReal;  
 
   Tu = (floor (time) - 0.5) / 36525;
   
   /* Greenwich sidereal time at 0 hrs UT (in radians) */
   sideRealG0 = 1.7399359 + 628.33195 * Tu + 0.00000675582 * square (Tu);

   sideRealG = sideRealG0 + (time - floor (time)) * SIDEREALSOLAR * TWOPI;

   localSideReal = sideRealG + longitude;

   return localSideReal;
}


