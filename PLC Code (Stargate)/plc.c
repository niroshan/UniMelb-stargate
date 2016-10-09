/*
 ********************************************************************************
 PLC.C

 Task to communicate to the PLC.

Author:	Len Sciacca
Date:	1994
 ********************************************************************************
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "plc.h"

/* #include "seq.h"   *//* Import definisitions, start_high, stop etc */


#define PLC_Fail 		0x01
#define PLC_Trip1		0x02
#define PLC_Trip		0x04
#define PLC_Spare_M99 	        0x08
#define PLC_24VDC   	        0x10
#define PLC_AzBrake		0x20
#define PLC_ElBrake		0x40
#define PLC_AzFinalCW 	        0x80

#define PLC_AzFinalCCW	        0x01
#define PLC_ElFinalUp	        0x02
#define PLC_ElFinalDwn	        0x04
#define PLC_EmergStop	        0x08
#define PLC_SmokeDetect	        0x10
#define PLC_Door1		0x20
#define PLC_Door2		0x40
#define PLC_CB5_6		0x80

#define PLC_CB4			0x01
#define PLC_CB3			0x02
#define PLC_RUN1		0x04
#define PLC_RUN2		0x08
#define PLC_C1			0x10
#define PLC_C2			0x20




/*
** PLC structure contains the actual fault/conditions
** PLC-bit contains the bit image of the PLC words
*/ 
static PLC_struct PLC_bit;


/*
 ****************************************************************************

 Data is encripted as follows
Word1: 	Bit 0 = Fail 	M96
Bit 1 = Trip1	        M97
Bit 2 = Trip	        M98
Bit 3 = Spare           M99
Bit 4 = 24VDC OK   	M100
Bit 5 = Az Brake	M101
Bit 6 = El Brake	M102
Bit 7 = Az Final CW 	M103

Word2:
Bit 0 = Az Final CCW	M104
Bit 1 = El Final UP	M105
Bit 2 = El Final DWN	M106
Bit 3 = Emerg. Stop	M107
Bit 4 = Smoke Detect	M108
Bit 5 = Door 1		M109
Bit 6 = Door 2		M110
Bit 7 = CB5/6		M111
word 3:
Bit 0 = CB4		M112
Bit 1 = CB3		M113
Bit 2 = RUN 1		M114
Bit 3 = RUN 2		M115
Bit 4 = C1		M116
Bit 5 = C2		M117
Bit 6 = Spare		M118
Bit 7 = Spare		M119

****************************************************************************
*/


int ReadPLC ( char word1, char word2, char word3, PLC_struct *PLC)
{
   PLC_bit.Fail      = (( PLC_Fail & word1 ) == PLC_Fail );
   PLC_bit.Trip1     = (( PLC_Trip1 & word1 ) == PLC_Trip1 );
   PLC_bit.Trip      = (( PLC_Trip & word1 ) == PLC_Trip );
   PLC_bit.C24VDC    = (( PLC_24VDC & word1 ) == PLC_24VDC );
   PLC_bit.AzBrake   = (( PLC_AzBrake & word1 ) == PLC_AzBrake );
   PLC_bit.ElBrake   = (( PLC_ElBrake & word1 ) == PLC_ElBrake );
   PLC_bit.AzFinalCW = (( PLC_AzFinalCW & word1 ) == PLC_AzFinalCW );

   PLC_bit.AzFinalCCW = (( PLC_AzFinalCCW & word2 ) == PLC_AzFinalCCW);
   PLC_bit.ElFinalUp  = (( PLC_ElFinalUp & word2 ) == PLC_ElFinalUp );
   PLC_bit.ElFinalDwn = (( PLC_ElFinalDwn & word2 ) == PLC_ElFinalDwn );
   PLC_bit.EmergStop  = (( PLC_EmergStop & word2 ) == PLC_EmergStop );
   PLC_bit.SmokeDetect	= (( PLC_SmokeDetect & word2 ) == PLC_SmokeDetect );
   PLC_bit.Door1	    = ((PLC_Door1 & word2 ) == PLC_Door1 );
   PLC_bit.Door2		= ((PLC_Door2 & word2 ) == PLC_Door2 );
   PLC_bit.CB5_6		= ((PLC_CB5_6 & word2 ) == PLC_CB5_6 );

   PLC_bit.CB4			= ((PLC_CB4 & word3 ) == PLC_CB4 );
   PLC_bit.CB3			= ((PLC_CB3 & word3 ) == PLC_CB3 );
   PLC_bit.RUN1		= ((PLC_RUN1 & word3 ) == PLC_RUN1 );
   PLC_bit.RUN2		= ((PLC_RUN2 & word3 ) == PLC_RUN2 );
   PLC_bit.C1			= ((PLC_C1 & word3 ) == PLC_C1 );
   PLC_bit.C2			= ((PLC_C2 & word3 ) == PLC_C2 );

   /*PLC = PLC_bit; */
   PLC->Fail      = PLC_bit.Fail;
   PLC->Trip1     = !PLC_bit.Trip1;
   PLC->Trip      = !PLC_bit.Trip;
   PLC->C24VDC    = !PLC_bit.C24VDC;
   PLC->AzBrake   = PLC_bit.AzBrake;
   PLC->ElBrake   = PLC_bit.ElBrake;
   PLC->AzFinalCW = !PLC_bit.AzFinalCW;

   PLC->AzFinalCCW = !PLC_bit.AzFinalCCW;
   PLC->ElFinalUp  = !PLC_bit.ElFinalUp;
   PLC->ElFinalDwn = !PLC_bit.ElFinalDwn;
   PLC->EmergStop  = !PLC_bit.EmergStop;
   PLC->SmokeDetect = PLC_bit.SmokeDetect;
   PLC->Door1	 = PLC_bit.Door1;
   PLC->Door2	 = PLC_bit.Door2;
   PLC->CB5_6	 = PLC_bit.CB5_6;

   PLC->CB4	 = !PLC_bit.CB4;
   PLC->CB3	 = !PLC_bit.CB3;
   PLC->RUN1	 = PLC_bit.RUN1;
   PLC->RUN2	 = PLC_bit.RUN2;
   PLC->C1		 = PLC_bit.C1;
   PLC->C2		 = PLC_bit.C2;

   
   return 1;
}

