CC = g++

COMMON_OBJECTS = obj/fft.o obj/misc.o obj/iadmm.o obj/admm1.o obj/admm05.o obj/matrix/cmat.o

DEBUG_COMMON_OBJECTS = debugObj/fft.o debugObj/misc.o debugObj/iadmm.o debugObj/admm1.o debugObj/admm05.o debugObj/matrix/cmat.o

# ALL_OBJECTS = $(COMMON_OBJECTS) obj/test.o obj/deblur.o
# DEBUG_ALL_OBJECTS = $(COMMON_OBJECTS) obj/test.o obj/deblur.o



LIBS = -lopencv_core -lopencv_highgui -lopencv_imgproc -lopencv_imgcodecs -lfftw3
CFLAGS = -O3 -std=c++11 -Wall
DEBUGFLAGS = -O0 -std=c++11 -g


debugObj/%.o: %.cpp
	$(CC) -c -o $@ $< $(DEBUGFLAGS)

%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)


test: $(COMMON_OBJECTS) obj/test.o
	$(CC) -o test $(COMMON_OBJECTS) obj/test.o $(LIBS)

debugTest: $(DEBUG_COMMON_OBJECTS) debugObj/test.o
	$(CC) -o debugTest $(DEBUG_COMMON_OBJECTS) debugObj/test.o $(LIBS)

deblur: $(COMMON_OBJECTS) obj/deblur.o
	$(CC) -o $@ $(COMMON_OBJECTS) obj/deblur.o $(LIBS)

%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

all:
	make test
	make deblur



.PHONY: clean

clean:
	rm -f $(ALL_OBJECTS) test deblur
