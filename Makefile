.PHONY: all clean

CFLAGS = -g -O2 -Wall -Wno-error -Wno-unused-variable -fsanitize=address
EXE= prototype4


all:
		gcc $(CFLAGS) -o $(EXE) cli.c blockjoin.c main.c -lz -L htslib/ -lhts -Wl,-rpath,htslib/ 
clean:
		rm -rf *.o $(EXE)
