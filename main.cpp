#include <iostream>
#include <opencv2/opencv.hpp>
#include <cmath>
#include <vector>
#include "matrix/cmat.hpp"
#include "fft.hpp"
#include "misc.hpp"



int main(int argc, char const *argv[]) {
  // RANDOM NUMBER GENERATOR FROM OPENCV LIB
  cv::RNG rng(12345);
  // char key;

  // std::cout << std::scientific;
  // std::cout << std::setprecision(10);

  int kernelSize = 17;
  double kernelSigma = 7;
  double sigma = 1.0e-6;
  double alpha=10;
  double mu = 0.05 / cv::max(sigma,1.e-12);




  cv::Mat cv_original;
  cv_original = cv::imread("../cameraman.tif", cv::IMREAD_GRAYSCALE);
  // cv::imshow("original", cv_original);
  // key = cv::waitKey(0);

  // LOADS ORIGINAL IMAGE TO A MAT
  Mat img(toMat(cv_original));

  // GENERATES KERNEL FOR BLURRING
  Mat k(kernel(kernelSize, kernelSigma));

  // GENERATES NOISE DATA
  Mat noiseData(img.rows, img.cols);
  fillNoise(noiseData, rng, 0, sigma);

  // BLUR image
  Mat blurred;
  convolute(blurred, img, k);

  // GENERATE NOISED MATRIX
  Mat noised;
  noised = img + noiseData;
  mapMat(noised, 0.0, 1.0);

  // GET CONSTANT MATRICES
  // CMat conjoDx, conjoDy, Nomin1, Denom1, Denom2;
  // getC (conjoDx, conjoDy, Nomin1, Denom1, Denom2, img, k);

  // START WITH THE CODE



  // % initialization
  // [m n]=size(Bn);
  // U=Bn; iter = 0;
  // beta = 1000;
  // Px=zeros(m,n);
  // Py=zeros(m,n);
  // PxB=Px;
  // PyB=Py;
  // S=zeros(100,1);

  // int rows = blurred.rows, cols = blurred.cols;
  // CMat U(blurred);
  // double gamma = 1000.0 / mu;
  // CMat Px(rows, cols);
  // CMat Py (rows, cols);
  // CMat PBx(rows, cols);
  // CMat PBy (rows, cols);
  //
  //
  // for(int k=0; k<10; ++k) {
  //
  // }




  // TESTING AREA


  // std::vector<Complex> v(5);
  // for(int i=0; i<5; ++i) v[i] = Complex(rand()%10, rand()%10);
  // std::cout << "v =";
  // for(int i=0; i<5; ++i)
  //   std::cout << "  " << v[i];
  // std::cout << std::endl;


  CMat a(3,3);
  CMat b(3,3);
  getRand(a);
  getRand(b);
  std::cout << "a = " << '\n' << a << std::endl;
  std::cout << "b = " << '\n' << b << std::endl;

  std::cout << "snr(a, b) = " << snr(a, b) << std::endl;
  // std::cout << "1/a = " << (a/a)/a << std::endl;





  return 0;
}
