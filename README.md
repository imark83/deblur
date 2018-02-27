## ADMM(1), ADMM(1/2) and IADMM

This repository contains the code in C++ to test codes to deblur images. In the deblurring process, the blurring operator is supposed to be known

This code includes a C++ wrapper for the library [fftw-3](http://www.fftw.org)
which this code deppends on. Also it uses [OpenCV Library](https://opencv.org) for loading and plotting images.


### Build commands
The following commands run on Unix machine. Be sure to have the [**FFTW-3 library**](http://www.fftw.org) and [**OpenCV library**](https://opencv.org) installed.

- To build the code for deblurring a single image:
```Makefile
make deblur
```
- To build the code to generate the test section data:
```Makefile
mkdir testData
make test
```

### Execute tests
Be sure you have the images **cameraman.tif**, **lena256.png** and **man256.png** in the same directory. (See *main.cpp:37* for details).
```
./tests
```
### Blur and deblur single images
```
./deblur [IMAGE NAME]
```
Type *./deblur* for help
