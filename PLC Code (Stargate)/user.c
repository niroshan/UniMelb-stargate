#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "general.h"
#include "user.h"


FILE *f_fromuser, *f_touser;

static double get_val (void);
static void get_name (char name[]);


int get_command (cp_data_struct *cp_data, unsigned int atcustate)
{
   char line[LINE];
    
   fflush (f_fromuser);

   fgets (line, LINE, f_fromuser);

   switch (line [0]) {
      case SATTRACK: {
         switch (atcustate) {
            case OFF_STATE:
               printf ("command is %c\n", line[0]);
               cp_data->cmd = POWER_ON_CMD;
            break;
            
            case STANDBY:
               printf ("command is %c\n", line[0]);
               cp_data->cmd = DRIVE_START_CMD;
            break;
            
            case HOLDING: {
               printf ("command is %c\n", line[0]);
               cp_data->cmd = GO_CMD;
               cp_data->mode = LEOTRACK_MODE;
               get_name (cp_data->name);
            }
            break;

            default:
               printf ("command is %c\n", line[0]);
            break;
         }
      }
      break;

      case STARTRACK: {
         switch (atcustate) {
            case OFF_STATE:
               printf ("command is %c\n", line[0]);
               cp_data->cmd = POWER_ON_CMD;
            break;
            
            case STANDBY:
               printf ("command is %c\n", line[0]);
               cp_data->cmd = DRIVE_START_CMD;
            break;
            
            case HOLDING: {
               printf ("command is %c\n", line[0]);
               cp_data->cmd = GO_CMD;
               cp_data->mode = STARTRACK_MODE;
               get_name (cp_data->name);
            }
            break;

            default:
               printf ("command is %c\n", line[0]);
            break;
         }
      }
      break;

      case MANUALTRACK: {
         switch (atcustate) {
            case OFF_STATE:
               printf ("command is %c\n", line[0]);
               cp_data->cmd = POWER_ON_CMD;
            break;
            
            case STANDBY:
               printf ("command is %c\n", line[0]);
               cp_data->cmd = DRIVE_START_CMD;
            break;
            
            case HOLDING: {
               printf ("command is %c\n", line[0]);
               cp_data->cmd = GO_CMD;
               cp_data->mode = MANUALTRACK_MODE;
               cp_data->az_cmd = get_val ();
               cp_data->el_cmd = get_val ();
            }
            break;

            default:
               printf ("command is %c\n", line[0]);
            break;
         }
      }
      break;
       
      case STOP: {
         if ((atcustate != OFF_STATE) 
               && (atcustate != STANDBY) 
               && (atcustate != STOPPING)) {
             printf ("command is %c\n", line[0]);
             cp_data->cmd = DRIVE_STOP_CMD;
         }

         else {
            printf ("command is %c\n", line[0]);
            cp_data->cmd = IDLE_CMD;
         }
      }
                 break;

      default:
         printf ("command is %c\n", line[0]);
      break;   
   }

   rewind (f_fromuser);

   return 1;
}


int send_feedback (double az, double el)
{
   int int_az, int_el;

   int_az = (int) az;
   int_el = (int) el;

   fprintf (f_touser, "*\n");
   fprintf (f_touser, "%d\n", int_az);
   fprintf (f_touser, ":\n");
   fprintf (f_touser, "%d\n", int_el);

   fflush (f_touser);

   return 1;
}


double get_val () 
{
   double val;

   fscanf (f_fromuser, "%lf", &val);

   return val;
}


void get_name (char name [])
{
   int i;

   fgets (name, NAME-1, f_fromuser);
   name [NAME-1] = NULLCHAR;

   for (i = 0; i < NAME; i++) {
      if (isspace(name[i]))
         break;
   }
   name[i] = NULLCHAR;
}
