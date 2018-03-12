CC = g++
COMMON_OBJECTS = fft.o misc.o iadmm.o admm1.o admm05.o lazyMat/cmat.o
ALL_OBJECTS = $(COMMON_OBJECTS) test.o deblur.o
LIBS = -lopencv_core -lopencv_highgui -lopencv_imgproc -lopencv_imgcodecs -lfftw3
CFLAGS = -O3 -std=c++11 -Wall

test: $(COMMON_OBJECTS) test.o
	$(CC) -o test $(COMMON_OBJECTS) test.o $(LIBS)

deblur: $(COMMON_OBJECTS) deblur.o
	$(CC) -o $@ $(COMMON_OBJECTS) deblur.o $(LIBS)

%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

all:
	make test
	make deblur



.PHONY: clean

clean:
	rm -f $(ALL_OBJECTS) test deblur
