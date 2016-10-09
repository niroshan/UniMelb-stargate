#CC=gnuc
CC=gcc
CFLAGS= -g -Wall -DDEBUG


OBJ=\
         error.o\
         iohw.o\
         main.o\
         packet.o\
         plc.o\
         poscon.o\
         sat_data.o\
         sat_track.o\
         seq.o\
         startrak.o\
         serial.o\
         user.o\
         utilities.o

all:tags prog

.PHONY : tags
tags: 
	ctags *.[ch]

prog: $(OBJ)
	$(CC) $(FLAGS) -o sequ $(OBJ) -lm

# TO back up the important filez
.PHONY : backup
backup:
	mount /a
	tar zcvf /a/backup.tar.gz *
	umount /a



.PHONY : clean
clean:
	/bin/rm -rf tags *.o core
clobber:
	/bin/rm -rf tags *.o  sequ core

# DO NOT DELETE
