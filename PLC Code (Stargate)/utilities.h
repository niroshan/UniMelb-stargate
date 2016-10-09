
/********************************************************************
**
** File: "utilities.h"
**
** Written: Thayaparan Thanabalasingham
**
** Last Modified: 12. 07. 98
**
********************************************************************/


long getDayNum (int year, int month, int day);
double getTime (void);
void getDate (long DayNum,int *Year, int *Month, int *Day);
void printTime (double t);

/* debugging */
void debug_going_down (void);
void debug_going_up (void);
void print_spaces (void);

/* dummies */
void enable (void);
void disable (void);
