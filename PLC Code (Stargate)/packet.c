#include <stdio.h>
#include <stdlib.h>

#include "packet.h"
#include "serial.h"
#include "general.h"

#define NULLCHAR '\0'



/* something doesnot work as supposed to.
** May be the protocol
** 29.07.98
*/
int pack (SendPacket packet)
{
   int i;
   unsigned char checksum = 0;
   unsigned char dacaz[2], dacel[2], poll;
   unsigned char send_buf[40] = "";

   /* these five might be the problem */
   dacaz[0] = (unsigned char) (packet.dacaz % 256);
   dacaz[1] = (unsigned char) (packet.dacaz / 256);
   
   dacel[0] = (unsigned char) (packet.dacel % 256);
   dacel[1] = (unsigned char) (packet.dacel / 256);

   poll = packet.pflag;

   send_buf[0] = STX;
   send_buf[1] = poll;
   send_buf[2] = dacaz[0];
   send_buf[3] = dacaz[1];
   send_buf[4] = dacel[0];
   send_buf[5] = dacel[1];
   send_buf[6] = packet.plc;
   send_buf[7] = NULLCHAR;

   for (i = 0; i < 7; i++) {
      checksum += (unsigned char) send_buf[i];
   }

   send_buf[7] = checksum;
   send_buf[8] = ETX;
   send_buf[9] = NULLCHAR;

   /*
   for (i = 0; i < 9; i++) {
      printf ("send_buf [i] = %d\n", send_buf [i]); 
   }
   */
   
   send_to_serial (send_buf, 9);

   return 0;    /* success */
}


int unpack (RcvdPacket *packet)
{
   int i;
   unsigned char azencoder[2], elencoder[2], adc[2];
   unsigned char word1, word2, word3;
   unsigned char checksum = 0;
   static unsigned char *rcvd = NULL;
   int length = 12;

   if (rcvd == NULL) {
      if ( (rcvd = (char *) malloc (80 * sizeof (char))) == NULL) {
          fprintf (stdout, "Cannot allocate memory\n"); 
          exit (1);
      }
   }
   /* ####### 
   **  Steps : 
   **     1. pack the packet in any form
   **     2. Call NCR's funciton which takes the packet and length
   **        & transmits it
   **
   */

   receive_from_serial (rcvd, length);  /* from NCR */

   /*
   printf ("Packet length = %d\n", strlen (rcvd));
   */
   if (length == 0) {
      printf ("Empty packet\n");
      return 2; /* no response */
   }
   
   for (i = 0; i < 10; i++) {
      checksum += (unsigned char) rcvd[i];
   }

   /*
   for (i=0;i<12;i++) {
      printf ("received: [%d]=%d\n", i,  rcvd[i]);
   }
   */

   if (rcvd[10] == checksum) {
      azencoder[0] = rcvd[1];
      azencoder[1] = rcvd[2];
      elencoder[0] = rcvd[3];
      elencoder[1] = rcvd[4];
      adc[0] = rcvd[5];
      adc[1] = rcvd[6];
      word1 = rcvd[7];
      word2 = rcvd[8];
      word3 = rcvd[9];
   }

   else {
      /*
      printf ("Check sum error *%d* vs *%d*\n", rcvd[10], checksum);
      */
      return 3; /* checksum error */
   }

   /* these might give problems : check the protocols */
   packet->azencoder = 256;
   packet->azencoder *= (int) azencoder[1];
   packet->azencoder += (int) azencoder[0];
   packet->elencoder = 1;
   packet->elencoder *= (int) azencoder[1];
   packet->elencoder += (int) azencoder[0];
   packet->adc = 1;
   packet->adc *= (int) adc[1];
   packet->adc += (int) adc[0];
   packet->word1 = word1;
   packet->word2 = word2;
   packet->word3 = word3;
   
   return 0;    /* success */
}

