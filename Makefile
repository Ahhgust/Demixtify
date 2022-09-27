
CC=g++
CFLAGS=-Wall -std=c++11 -Ofast -DDEBUG=0 -fomit-frame-pointer
#CFLAGS=-Wall -std=c++11 -ggdb3 -DDEBUG=1
#CFLAGS=-Wall -std=c++11 -ggdb3

# taken from stack exchange!
SYS := $(shell ${CC} -dumpmachine)


# TODO:
#htslib must be enabled w/o libcurl (cannot be combined with -static) :
# ./configure --disable-libcurl


ifneq (, $(findstring mingw, $(SYS)))
# note2self: -static only works if disable-libcurl is flagged w/ libcurl
# also... pthreads? are we sure?
ALL: demix.o
	${CC} ${CFLAGS} demix.o -o demix.exe  htslib/libhts.a  -lz -lbz2 -llzma -lm -static
else

ALL: demix.o
	${CC} ${CFLAGS} demix.o -o demix  htslib/libhts.a  -lz -lbz2 -llzma -lpthread -lm -ldeflate
endif


demix.o: demix.c demix.h
	${CC} ${CFLAGS} -c demix.c


ifneq (, $(findstring mingw, $(SYS)))
clean:
	del demix.exe demix.o
else
clean:
	rm demix demix.o
endif
