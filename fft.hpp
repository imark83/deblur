#ifndef __FFT_HPP__
#define __FFT_HPP__

#include <complex>
#include <vector>
#include <fftw3.h>

#include "matrix/cmat.hpp"
typedef std::vector<std::complex<double> > CArray;


void fft(CArray &x);
void ifft(CArray &x);
void fft(CMat &mat);
void ifft(CMat &mat);
void fftn(CMat &mat);
void ifftn(CMat &mat);

#endif
