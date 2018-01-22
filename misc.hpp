#ifndef __MISC_HPP__
#define __MISC_HPP__

#include <opencv2/opencv.hpp>
#include "matrix/cmat.hpp"


// FILLS WITH GAUSSIAN NOISE
void fillNoise(Mat &rop, cv::RNG &rng, double mean, double std);

// CUTS VALUES OF MATRIX TO MINVAL AND MAXMAL
void mapMat(Mat &rop, double minval, double maxval);


// CONVOLUTES MATRIX APPLYING KERNEL USING CIRCULAR MAPPING(TORUS)
void convolute(Mat &rop, const Mat &op, const Mat &kernel);

// PS2PTF FUNCTION
void psf2otf (CMat &rop, const CMat &op1, int outsize[2]);

// INITIALIZES CONSTANT STUFF
void getC(CMat &conjoDx, CMat &conjoDy, CMat &Nomin1, CMat &Denom1,
      CMat &Denom2, const Mat &Bn, const Mat &H);

// SHRINFT FUNCTION
void shrinft(CMat &x, double nameta);

// GENERATES GAUSSIAN KERNEL
// size MUST be odd
Mat kernel(int size, double sigma);

// CONVERTS CVMATRIX TO MATRIX
Mat toMat(const cv::Mat &op);
cv::Mat toCVMat(const Mat &op);
Mat toMat(const cv::Mat &op, double);
cv::Mat toCVMat(const Mat &op, double);



// MISC STATISTIC FUNCTIONS
template <class T>
T mean(T *op, int n) {
  T rop = T(0);
  for(int i=0; i<n; ++i)
    rop = rop + op[i];
  return (rop / ((double) n));
}

template <class T>
T mean(const Mat_<T> &op) {
  return mean<T>(op.data, op.rows*op.cols);
}

template <class T>
T var(T *op, int n) {
  T rop = T(0);
  T m = mean<T>(op, n);
  for(int i=0; i<n; ++i)
    rop = rop + (op[i] * op[i]);
  return ((rop / ((double) n)) - m*m);
}

template <>
Complex var(Complex *op, int n);


template <class T>
T var(const Mat_<T> &op) {
  return var(op.data, op.rows*op.cols);
}



// FUNCTION SNR
template <class T>
T snr(const Mat_<T> &ref, const Mat_<T> &sig) {
  T mse = (mean(pow2(ref-sig)));
  T dv = var(ref);
  return 10.0 * log10(dv/mse);
}

/*function x = snr(ref, sig)
% snr -- Compute Signal-to-Noise Ratio for images
%
% Usage:
%       x = snr(ref, sig)
%
% Input:
%       ref         Reference image
%       sig         Modified image
%
% Output:
%       x           SNR value

mse = mean((ref(:)-sig(:)).^2);
dv = var(ref(:),1);
x = 10*log10(dv/mse);
*/


#endif
