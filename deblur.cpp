#include <iostream>
#include <fstream>
#include <iomanip>
#include <opencv2/opencv.hpp>
#include <cmath>
#include <vector>
#include "matrix/cmat.hpp"
#include "fft.hpp"
#include "misc.hpp"
#include "admm1.hpp"
#include "admm05.hpp"
#include "iadmm.hpp"

#include <stdio.h>

void toTXT (Mat &op, const char *fname) {
  std::ofstream fout;
  fout.open(fname);
  fout << std::scientific << std::setprecision(16);
  for(int i=0; i<op.rows; ++i) {
    fout << op(i,0);
    for(int j=1; j<op.cols; ++j)
      fout << " " << op(i,j);
    fout << std::endl;
  }
  fout.close();
  return;
}


void initParameters (int &kernelSize, double &kernelSigma,
      double &sigma, double &alpha, int &nIter,
      int argc, const char *argv[]) {

  if (argc > 2)
    nIter = atoi(argv[2]);
  if (argc > 3)
    kernelSize = atoi(argv[3]);
  if (argc > 4)
    kernelSigma = atof(argv[4]);
  if (argc > 5)
    sigma = atof(argv[5]);
  if (argc > 6)
    alpha = atof(argv[6]);
  return;
}


void help(const char *argv[]) {
  printf("Ussage:\n%s IMAGE_NAME [nIter (200)]"
          "[kernelSize (17)] [kernelSigma (7)]"
          " [sigma (1.0e-6)]" " [alpha (10)] \n", argv[0]);
  exit(1);
}

int main(int argc, char const *argv[]) {
  // RANDOM NUMBER GENERATOR FROM OPENCV LIB
  cv::RNG rng(12345);

  std::cout << std::scientific;
  std::cout << std::setprecision(10);

  if (argc == 1)
    help(argv);

  int kernelSize = 17;
  double kernelSigma = 7;
  double sigma = 1.0e-5;
  double alpha = 0.1;
  int nIter = 200;
  double mu = 1.0e9; //0.05 / cv::max(sigma,1.e-12);
  initParameters (kernelSize, kernelSigma, sigma, alpha, nIter, argc, argv);

  std::cout << "nIter = " << nIter << std::endl;
  std::cout << "kernelSize = " << kernelSize << std::endl;
  std::cout << "kernelSigma = " << kernelSigma << std::endl;
  std::cout << "sigma = " << sigma << std::endl;
  std::cout << "alpha = " << alpha << std::endl;

  cv::Mat cv_original;
  cv_original = cv::imread(argv[1], cv::IMREAD_GRAYSCALE);

#ifndef DONT_BLUR_ME
  cv::imshow("Original", cv_original);
  cv::moveWindow("Original", 50, 50);
#endif

  // LOADS ORIGINAL IMAGE TO A MAT
  Mat img(toMat(cv_original));

  // GENERATES KERNEL FOR BLURRING
  Mat k(kernel(kernelSize, kernelSigma));

  // GENERATES NOISE DATA
  Mat noiseData(img.rows, img.cols);
  fillNoise(noiseData, rng, 0, sigma);

  // BLUR image
  Mat blurred(img);
#ifndef DONT_BLUR_ME
  convolute(blurred, img, k);
#endif

  // GENERATE NOISED MATRIX
  Mat noised(blurred);
#ifndef DONT_BLUR_ME
  noised = blurred + noiseData;
  mapMat(noised, 0.0, 1.0);
#endif

  cv::imshow("Blurred", toCVMat(noised));
  cv::moveWindow("Blurred", img.cols + 100, 50);
  cv::waitKey(100);


  Mat rop;
  std::vector<double> OBJ, E, TV, S;

  rop = admm1(OBJ, TV, E, S, img, k, noised, mu, alpha, nIter);

  cv::imshow("Deblurred ADMM(1)", toCVMat(rop));
  cv::moveWindow("Deblurred ADMM(1)", 2*img.cols + 150, 50);
  cv::waitKey(100);



  rop = admm05(OBJ, TV, E, S, img, k, noised, mu, alpha, nIter);

  cv::imshow("Deblurred ADMM(1/2)", toCVMat(rop));
  cv::moveWindow("Deblurred ADMM(1/2)", 3*img.cols + 200, 50);
  cv::waitKey(100);


  rop = iadmm(OBJ, TV, E, S, img, k, noised, mu, alpha, nIter);

  cv::imshow("Deblurred IADMM", toCVMat(rop));
  cv::moveWindow("Deblurred IADMM", 4*img.cols + 250, 50);
  cv::waitKey(0);

  // TESTING AREA


  // std::vector<Complex> v(5);
  // for(int i=0; i<5; ++i) v[i] = Complex(rand()%10, rand()%10);
  // std::cout << "v =";
  // for(int i=0; i<5; ++i)
  //   std::cout << "  " << v[i];
  // std::cout << std::endl;


  // CMat a(3,3);
  // CMat b(3,3);
  // getRand(a);
  // getRand(b);
  // std::cout << "a = " << '\n' << a << std::endl;
  // std::cout << "b = " << '\n' << b << std::endl;
  //
  // std::cout << "snr(a, b) = " << snr(a, b) << std::endl;
  // std::cout << "1/a = " << (a/a)/a << std::endl;





  return 0;
}
