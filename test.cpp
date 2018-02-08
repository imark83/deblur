#include <iostream>
#include <fstream>

#include <opencv2/opencv.hpp>
#include "matrix/cmat.hpp"
#include "misc.hpp"
#include "iadmm.hpp"

std::ostream &operator<<(std::ostream &output, const CArray &op) {
  for(size_t i=0; i<op.size(); ++i)
    output << i+1 << "  " << op[i].real() << std::endl;
  return output;
}
std::ostream &operator<<(std::ostream &output, const std::vector<double> &op) {
  for(size_t i=0; i<op.size(); ++i)
    output << i+1 << "  " << op[i] << std::endl;
  return output;
}

int main(int argc, char const *argv[]) {
  // RANDOM NUMBER GENERATOR FROM OPENCV LIB
  cv::RNG rng(1);

  // PARAMETERS
  int kernelSize = 17;
  double kernelSigma = 7;
  double sigma = 1.0e-6;
  double alpha=10;
  int nIter = 1000;
  double mu = 0.05 / cv::max(sigma,1.e-12);

  cv::Mat cv_original1 = cv::imread("../cameraman.tif", cv::IMREAD_GRAYSCALE);
  // LOADS ORIGINAL IMAGE TO A MAT
  Mat I1(toMat(cv_original1));

  // GENERATES KERNEL FOR BLURRING
  Mat H(kernel(kernelSize, kernelSigma));

  // BLUR IMAGE WITH GAUSSIAN FILTER
  Mat B1;
  convolute(B1,I1,H);

  // GENERATES NOISE DATA
  Mat noiseData(I1.rows, I1.cols);
  fillNoise(noiseData, rng, 0, sigma);

  // GENERATE NOISED MATRIX
  Mat Bn1(I1.rows, I1.cols);
  Bn1 = B1 + noiseData;
  mapMat(Bn1, 0.0, 1.0);


  Mat rop;
  CArray S1;
  std::vector<double> E1;
  rop = iadmm(E1, S1, I1, H, Bn1, mu, alpha, nIter);


  std::ofstream fout;
  fout.open("testData/s1.txt");
  fout << S1 << std::endl;
  fout.close();

  fout.open("testData/e1.txt");
  fout << E1 << std::endl;
  fout.close();


  return 0;
}
