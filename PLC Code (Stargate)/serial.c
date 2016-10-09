#include <stdio.h>
#include <stdlib.h>
#include <sys/syscall.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <sys/ioctl.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <termios.h>
#include <errno.h>

#define COM_PORT  "/dev/ttyS0" 
/*
#define COM_PORT  "/dev/ttyS3" 
*/

#define REPETITION      3

#define TRUE    1
#define FALSE   0
#define SER_DATA 1      /* DATA packet        */


int serial_handle;
struct termios oldstty;


static void Setuptty(int file, struct termios *currterm);
static int is_identical (const char *str, size_t start,  size_t num);
static void print_packet (const char *pak);
static int is_valid (const char *msg);
static void process_packet (char *pak);


/* Sets up the terminal associated with the file descriptor file,
 ** and returns the values required to reset it in *currterm.
 */
void Setuptty(int file, struct termios *currterm)
{
   char tbuf[500];
   struct termios newterm;

   ioctl(file, TCGETS, currterm);
   newterm=(*currterm);

   newterm.c_oflag &= (~OPOST);
   newterm.c_lflag = newterm.c_lflag & (~ICANON);

   
   if (ioctl(file, TCSETSW, &newterm) != 0) {
      printf ("ioctl error \n");        /* Thaya */
      exit (1);
   }
   
   /* This is very kludgy, should (nah, WILL) fix so it uses the ioctl above */
   sprintf(tbuf, "stty 0:4:cbd:0:3:1c:7f:15:4:0:1:0:11:13:1a:0:12:f:17:16:0:0:90 < %s", COM_PORT);
   system(tbuf);

   return;
}


void init_serial ()
{
   if ((serial_handle=open(COM_PORT,O_RDWR | O_NDELAY ))== -1) {

      printf ("Error: Couldn't open serial line!\n");
      exit (1);
   }
   Setuptty(serial_handle,&oldstty); /* Set up the serial port */      
   tcflush(serial_handle, TCIFLUSH);  
}

/* send_to_serial(buffer,datasize)
 * Tries to write datasize bytes to the serial line, starting
 * at buffer. If for some reason the bytes may not be written,
 * the function keeps trying until the maximum amount of
 * time has elapsed.
 */
int send_to_serial(unsigned char *buffer, int datasize) {
   int i,wrote,numfds;
   fd_set writefds;
   struct timeval writeto;

   /* Thaya 13. 08. 98 */
   printf ("Message sent: %s\n", buffer);

   tcflush (serial_handle, TCOFLUSH);


   writeto.tv_sec=5;
   writeto.tv_usec=0;
   i=wrote=0;
   numfds=1;

   FD_ZERO(&writeto);
   FD_SET(serial_handle,&writeto);

   while((i < datasize)&&(numfds>0)) {
      wrote=write(serial_handle,buffer + i,datasize - i);
      i+=wrote;
      if (i<datasize)
         numfds=select(serial_handle+1,(fd_set *)NULL,
               &writefds, (fd_set *)NULL, &writeto);
   }

   return 1;
}


/* receive_from_serial(buffer,datasize)
 * Tries to read datasize bytes of data from the serial
 * connection, starting at buffer. If no data is available
 * it will try to read for a given amount of time. If no
 * data is still available it will return -1. If data was
 * read successfully it will return 1.
 */
int receive_from_serial(char *buffer, int datasize) {
   int i,red,numfds;
   fd_set readfds;
   struct timeval readto;

   readto.tv_sec=4;
   readto.tv_usec=0;
   i=red=0;

   FD_ZERO(&readfds);
   FD_SET(serial_handle,&readfds);
   numfds=1;
   printf ("Just Before reading\n");
   
   while((i<datasize)&&(numfds>0)) {
      red=read(serial_handle,buffer + i,datasize - i);
      printf ("Just after read () call.\n");
      if (red>0)
         i+=red;
      if (i<datasize)
         numfds = select(serial_handle+1,&readfds,
               (fd_set *)NULL, (fd_set *)NULL, &readto);
   }
   buffer[i] = '\0';

   if (numfds==0) {
      printf ("Timed out for a packet\n");
      numfds=(-1);
   }
   /* Thaya 13. 08. 98*/
   else {
      printf ("Message: %s\n", buffer);
   }
    
   tcflush (serial_handle, TCIFLUSH);

   /* just now - thaya 13. 08. 98
   process_packet (buffer);
   */
   return numfds;
}


/**  is_identical ():
 *      checks if <NUM> subsequent chars from <START> are the same
 *
 *   returns TRUE/ FALSE
 */
int is_identical (const char *str, size_t start,  size_t num)
{
   int i;
   for (i=0 ; i < num-1 ; i++)
      if (str[i+start] != str [i+start+1])
         return FALSE;

   return TRUE;

}


void print_packet (const char *msg)
{
   int i;

   fprintf (stderr, "Printing packet!\n");
   for (i=0 ; i < strlen (msg) ; i++)
      fprintf (stderr, "\tmsg[%d] = %d\n", i, (int) msg[i]);
   fprintf (stderr, "\n");


}


/** is_valid ():
 ** ################## need to more error checking ################
 **     checks if the msg is intact
 **
 **  returns: TRUE/FALSE
 **/
int is_valid (const char *msg)
{
   size_t packet_len = 0;
   size_t packet_num = 0;

   /*****  process data packet   */

   /* 1. First char is SER_DATA  */
   size_t curr_loc = 0;

   if ( is_identical (msg, 0, REPETITION) == FALSE)
      return FALSE;
   curr_loc += REPETITION;

   /* 2. packet number */
   packet_num = msg [curr_loc];
   curr_loc += REPETITION;

   /* 3. ###### Pack length 
    *  STill need to be decoded.  Assuming 1byte (0-255) now
    */
   packet_len = (size_t) msg [curr_loc];
   curr_loc += REPETITION;


   return TRUE;

}


void process_packet (char *pak)
{
   char *start;

   if (is_valid (pak) == FALSE)
   {
      fprintf (stderr, "packet not valid!.. Dropping it!!\n");
      return;
   }

   /* by pass all the header info in pac */
   start = (pak + 3 * REPETITION);

   /* strip off the last 2 CRC bytes */
   pak [strlen (start) - 2 ] = '\0';

}

