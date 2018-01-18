CC = g++
TARGETS = main.cpp fft.cpp misc.cpp matrix/cmat.cpp
OBJECTS = main.o fft.o misc.o matrix/cmat.o
LIBS = -lopencv_core -lopencv_highgui -lopencv_imgproc -lopencv_imgcodecs -lfftw3
CFLAGS = -std=c++11 -g -Wall


%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

aaa : $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LIBS)

all:
	make aaa



.PHONY: clean

clean:
	rm -f $(OBJECTS)
