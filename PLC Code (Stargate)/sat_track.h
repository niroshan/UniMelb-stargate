
/*********************************************************
**
** File: "sat_track.h"
**
** Written: Thayaparan Thanabalasingham
**
** Last Modified: 29. 07. 98 - 21:04:00
**
*********************************************************/

#include "seq.h"

typedef double mat3x3[3][3];


void satTrack (const char *s, double currentTime,
      const SiteStruct *site, double *azimuth, double *elevation);
