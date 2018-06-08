CC = g++

COMMON_OBJECTS = obj/fft.o obj/misc.o obj/iadmm.o obj/admm1.o obj/admm05.o obj/matrix/cmat.o

DEBUG_COMMON_OBJECTS = obj/debug/fft.o obj/debug/misc.o obj/debug/iadmm.o obj/debug/admm1.o obj/debug/admm05.o obj/debug/matrix/cmat.o

ALL_OBJECTS = $(COMMON_OBJECTS) obj/test.o obj/deblur.o
DEBUG_ALL_OBJECTS = $(DEBUG_COMMON_OBJECTS) obj/debug/test.o obj/debug/deblur.o



LIBS = -lopencv_core -lopencv_highgui -lopencv_imgproc -lopencv_imgcodecs -lfftw3
CFLAGS = -O3 -std=c++11 -Wall
DEBUGFLAGS = -O0 -std=c++11 -g

test: $(COMMON_OBJECTS) obj/test.o
	$(CC) -o test $(COMMON_OBJECTS) obj/test.o $(LIBS)

obj/debug/%.o: %.cpp
	$(CC) -c -o $@ $< $(DEBUGFLAGS)

obj/%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)


debugTest: $(DEBUG_COMMON_OBJECTS) obj/debug/test.o
	$(CC) -o debugTest $(DEBUG_COMMON_OBJECTS) obj/debug/test.o $(LIBS)

deblur: $(COMMON_OBJECTS) obj/deblur.o
	$(CC) -o $@ $(COMMON_OBJECTS) obj/deblur.o $(LIBS)

%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

all: clean
	make test
	make debugTest



.PHONY: clean

clean:
	rm -f $(ALL_OBJECTS) $(DEBUG_ALL_OBJECTS) test deblur debugTest
