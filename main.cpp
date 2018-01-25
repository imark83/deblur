#include <iostream>
#include <fstream>
#include <iomanip>
#include <opencv2/opencv.hpp>
#include <cmath>
#include <vector>
#include "matrix/cmat.hpp"
#include "fft.hpp"
#include "misc.hpp"
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
    kernelSize = atoi(argv[2]);
  if (argc > 3)
    kernelSigma = atof(argv[3]);
  if (argc > 4)
    sigma = atof(argv[4]);
  if (argc > 5)
    alpha = atof(argv[5]);
  if (argc > 6)
    nIter = atoi(argv[6]);
  return;
}


void help(const char *argv[]) {
  printf("Ussage:\n%s IMAGE_NAME [kernelSize (17)]"
          " [kernelSigma (7)] [sigma (1.0e-6)]"
          " [alpha (10)] [nIter (200)]\n", argv[0]);
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
  double sigma = 1.0e-6;
  double alpha=10;
  int nIter = 200;
  double mu = 0.05 / cv::max(sigma,1.e-12);
  initParameters (kernelSize, kernelSigma, sigma, alpha, nIter, argc, argv);

  std::cout << "kernelSize = " << kernelSize << std::endl;
  std::cout << "kernelSigma = " << kernelSigma << std::endl;
  std::cout << "sigma = " << sigma << std::endl;
  std::cout << "alpha = " << alpha << std::endl;
  std::cout << "nIter = " << nIter << std::endl;

  cv::Mat cv_original;
  cv_original = cv::imread(argv[1], cv::IMREAD_GRAYSCALE);
  cv::imshow("Original", cv_original);
  cv::moveWindow("Original", 50, 50);

  // LOADS ORIGINAL IMAGE TO A MAT
  Mat img(toMat(cv_original));

  // GENERATES KERNEL FOR BLURRING
  Mat k(kernel(kernelSize, kernelSigma));

  // GENERATES NOISE DATA
  Mat noiseData(img.rows, img.cols);
  fillNoise(noiseData, rng, 0, sigma);

  // BLUR image
  Mat blurred(img);
  convolute(blurred, img, k);
  cv::imshow("Blurred", toCVMat(blurred));
  cv::moveWindow("Blurred", 350, 50);
  cv::waitKey(100);

  // GENERATE NOISED MATRIX
  Mat noised;
  noised = img + noiseData;
  mapMat(noised, 0.0, 1.0);


  Mat rop;
  CArray S;

  rop = iadmm(S, img, k, blurred, mu, alpha, nIter);


  cv::imshow("Deblurred", toCVMat(rop));
  cv::moveWindow("Deblurred", 650, 50);
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
