#ifndef _GENERAL_H
#define _GENERAL_H
/************************************************************************/
/*									*/
/*				GENERAL.H				*/
/*									*/
/*		Collection of some commonly used definitions		*/
/*			for use by all modules				*/
/*									*/
/*	Programmer:	P. Moylan					*/
/*	Function:	This file contains definitions which		*/
/*			are needed in almost all modules.		*/
/*									*/
/************************************************************************/

#define TRUE    1
#define FALSE   0

#define STX     0x02
#define ETX     0x03

#define NULLCHAR        '\0'

#define LINE    80

#ifndef PI
#define PI  3.141592653589793238462643
#endif
#define DEGREES_PER_RADIAN (180 / PI)
#define RADIANS_PER_DEGREE (PI / 180)

/* for debugging */
extern int debug_spaces;

#endif /* __GENERAL_H */
