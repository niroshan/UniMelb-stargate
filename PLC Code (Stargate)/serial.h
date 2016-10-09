/*
** "02. 08.98"
** This contains the functions that are to be supplied by the ncr's
** module.
*/


void init_serial (void);
int send_to_serial (const unsigned char *send_buf, int len);
int receive_from_serial (unsigned char *buf, int len);
