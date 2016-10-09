typedef struct packet_to_send_struct {
   int dacaz;
   int dacel;
   unsigned char plc;
   int pflag;
} SendPacket;


typedef struct packet_rcvd_struct {
   int azencoder;
   int elencoder;
   int adc;
   unsigned char word1;
   unsigned char word2;
   unsigned char word3;
} RcvdPacket;


int pack (SendPacket packet);
int unpack2 (SendPacket *packet);
int unpack (RcvdPacket *packet);
int send (const char send_buf [], int len);
int receive (char buf [], int *len);
