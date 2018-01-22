#include "fft.hpp"


// WRAPPER FOR FFTW3 LIBRARY
void fft (std::vector<Complex> &vec) {
  fftw_complex *in;
  // fftw_complex *out;
  fftw_plan p;


  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * vec.size());
  for(size_t i=0; i<vec.size(); ++i) {
    in[i][0] = vec[i].real();
    in[i][1] = vec[i].imag();
  }
  // out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * vec.size());
  p = fftw_plan_dft_1d(vec.size(), in, in, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p); /* repeat as needed */

  for(size_t i=0; i<vec.size(); ++i) {
    vec[i] = Complex(in[i][0], in[i][1]);
  }

  fftw_destroy_plan(p);
  fftw_free(in);
  // fftw_free(out);

  return;
}


void ifft (std::vector<Complex> &vec) {
  fftw_complex *in;
  // fftw_complex *out;
  fftw_plan p;


  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * vec.size());
  for(size_t i=0; i<vec.size(); ++i) {
    in[i][0] = vec[i].real();
    in[i][1] = vec[i].imag();
  }
  // out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * vec.size());
  p = fftw_plan_dft_1d(vec.size(), in, in, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p); /* repeat as needed */

  for(size_t i=0; i<vec.size(); ++i) {
    vec[i] = Complex(in[i][0], in[i][1]) / ((double) vec.size());
  }

  fftw_destroy_plan(p);
  fftw_free(in);
  // fftw_free(out);

  return;
}

void fft(CMat &mat) {
  CArray x(mat.rows);
  for (int j=0; j<mat.cols; ++j) {
    for (int i=0; i<mat.rows; ++i) x[i] = mat(i,j);
    fft(x);
    for (int i=0; i<mat.cols; ++i) mat(i,j) = x[i];
  }
}

void ifft(CMat &mat) {
  CArray x(mat.rows);
  for (int j=0; j<mat.cols; ++j) {
    for (int i=0; i<mat.rows; ++i) x[i] = mat(i,j);
    ifft(x);
    for (int i=0; i<mat.cols; ++i) mat(i,j) = x[i];
  }
}

void fftn(CMat &mat) {
  fft(mat);
  std::complex<double> aux;
  for (int i=0; i<mat.rows; ++i) for (int j=i+1; j<mat.cols; ++j) {
    aux = mat(i,j); mat (i,j) = mat(j,i); mat(j,i) = aux;
  }
  fft(mat);
  for (int i=0; i<mat.rows; ++i) for (int j=i+1; j<mat.cols; ++j) {
    aux = mat(i,j); mat (i,j) = mat(j,i); mat(j,i) = aux;
  }
}

void ifftn(CMat &mat) {
  ifft(mat);
  std::complex<double> aux;
  for (int i=0; i<mat.rows; ++i) for (int j=i+1; j<mat.cols; ++j) {
    aux = mat(i,j); mat (i,j) = mat(j,i); mat(j,i) = aux;
  }
  ifft(mat);
  for (int i=0; i<mat.rows; ++i) for (int j=i+1; j<mat.cols; ++j) {
    aux = mat(i,j); mat (i,j) = mat(j,i); mat(j,i) = aux;
  }
}
