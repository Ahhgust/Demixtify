
CC=g++
#CFLAGS=-Wall -std=c++11 -Ofast -DDEBUG=0 -fomit-frame-pointer -DC11THREADS

CFLAGS=-Wall -std=c++11 -Ofast -DDEBUG=0 -fomit-frame-pointer

#CFLAGS=-Wall -std=c++11 -ggdb3 -DDEBUG=1

#CFLAGS=-Wall -std=c++11 -pg -ggdb3

# taken from stack exchange!
SYS := $(shell ${CC} -dumpmachine)

# Here is the static compile command.
# g++ -Wall -std=c++11 -Ofast -DDEBUG=0 -fomit-frame-pointer demix.o -o demix_static -L/usr/lib64 htslib/libhts.a  -lz -lbz2 -llzma -lpthread -lm -ldeflate -static
# Note:
# you need static libraries of everything. Centos/rocky linux makes that annoyingly hard.
# See: https://forums.centos.org/viewtopic.php?t=52129


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
