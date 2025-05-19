.PHONY: all clean

CFLAGS = -g -O2 -Wall -Wno-error -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-label 
EXE=pomfret
DIR_HTSLIB?=/usr/bin/htslib


all:
		gcc $(CFLAGS) -o $(EXE) kthread.c kstring.c cli.c blockjoin.c main.c -lz -pthread -L $(DIR_HTSLIB) -lhts -Wl,-rpath='$(DIR_HTSLIB)' -I $(DIR_HTSLIB)

dbg: CFLAGS+=-fsanitize=address
dbg: all

testbam:
		gcc $(CFLAGS) -o test test_modify_bam.c cli.c -fsanitize=address -L htslib/ -lhts -Wl,-rpath,$(DIR_HTSLIB)  -I $(DIR_HTSLIB)

clean:
		rm -rf *.o $(EXE)
