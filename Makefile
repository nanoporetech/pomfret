UNAME := $(shell uname -s)

.PHONY: all clean

CFLAGS = -g -O2 -Wall -Wno-error -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-label 
EXE=pomfret
DIR_HTSLIB?=/usr/bin/htslib


all:
		gcc $(CFLAGS) -o $(EXE) kthread.c kstring.c cli.c blockjoin.c main.c -lz -pthread -L $(DIR_HTSLIB) -lhts -Wl,-rpath,'$(DIR_HTSLIB)' -I $(DIR_HTSLIB)
ifeq ($(UNAME), Darwin)  # macos
	install_name_tool -id @rpath/libhts.3.dylib $(DIR_HTSLIB)/libhts.dylib
	install_name_tool -change /usr/local/lib/libhts.3.dylib libhts.dylib pomfret
endif

dbg: CFLAGS+=-fsanitize=address
dbg: all

clean:
		rm -rf *.o $(EXE)
