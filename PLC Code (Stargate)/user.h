#include "seq.h"

#define SATTRACK        'S'
#define STARTRACK       'R'
#define MANUALTRACK     'M'
#define STOP            'T'

extern FILE *f_fromuser, *f_touser;


int get_command (cp_data_struct *cpdata, unsigned int atcustate);
int send_feedback (double az, double el);
