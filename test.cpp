#include <iostream>
#include <fstream>
#include <cstdio>

#include <opencv2/opencv.hpp>
#include "matrix/cmat.hpp"
#include "misc.hpp"
#include "admm1.hpp"
#include "admm05.hpp"
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
  int nImages = 3;
  for(int k=0; k<nImages; ++k) {
    // RANDOM NUMBER GENERATOR FROM OPENCV LIB
    cv::RNG rng(1);

    // PARAMETERS
    int kernelSize = 17;
    double kernelSigma = 7;
    double sigma = 1.0e-6;
    double alpha=0.1;
    int nIter = 100;
    double mu = 1e12; //0.05 / cv::max(sigma,1.e-12);

    char fname[100];
    string imageName[] = {"cameraman.tif", "lena256.png", "man256.png"};

    cv::Mat cv_original = cv::imread(imageName[k], cv::IMREAD_GRAYSCALE);
    // LOADS ORIGINAL IMAGE TO A MAT
    Mat I(toMat(cv_original));

    cv::imshow("Original", cv_original);
    cv::moveWindow("Original", 50, 50);

    // GENERATES KERNEL FOR BLURRING
    Mat H(kernel(kernelSize, kernelSigma));

    // BLUR IMAGE WITH GAUSSIAN FILTER
    Mat B;
    convolute(B,I,H);

    // GENERATES NOISE DATA
    Mat noiseData(I.rows, I.cols);
    fillNoise(noiseData, rng, 0, sigma);

    // GENERATE NOISED MATRIX
    Mat Bn(I.rows, I.cols);
    Bn = B + noiseData;
    mapMat(Bn, 0.0, 1.0);

    cv::imshow("Blurred", toCVMat(Bn));
    cv::moveWindow("Blurred", I.cols + 100, 50);
    cv::waitKey(100);


    Mat rop;
    std::vector<double> OBJ, E, TV, S;
    std::ofstream fout;

    rop = admm1(OBJ, TV, E, S, I, H, Bn, mu, alpha, nIter);
    sprintf(fname, "testData/admm(1)-%02i.txt", k);
    fout.open(fname);
    for(int i=0; i<nIter; ++i)
      fout << i << " " << OBJ[i] << " " << TV[i]
            << " " << E[i] << " " << S[i] << std::endl;
    fout.close();

    rop = admm05(OBJ, TV, E, S, I, H, Bn, mu, alpha, nIter);
    sprintf(fname, "testData/admm(1/2)-%02i.txt", k);
    fout.open(fname);
    for(int i=0; i<nIter; ++i)
      fout << i << " " << OBJ[i] << " " << TV[i]
            << " " << E[i] << " " << S[i] << std::endl;
    fout.close();


    rop = iadmm(OBJ, TV, E, S, I, H, Bn, mu, alpha, nIter);

    sprintf(fname, "testData/iadmm-%02i.txt", k);
    fout.open(fname);
    for(int i=0; i<nIter; ++i)
      fout << i << " " << OBJ[i] << " " << TV[i]
            << " " << E[i] << " " << S[i] << std::endl;
    fout.close();

  }

  return 0;
}
