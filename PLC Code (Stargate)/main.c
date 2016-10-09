#include <stdio.h>
#include <unistd.h>

#include "seq.h"

int debug_spaces = 0;

int main (void)
{

   init_sequencer_task ();

   while (1) {
      sequencer_task ();
      usleep (40000); /* 100,000 us = 100ms */
   }
   return 0;
}
