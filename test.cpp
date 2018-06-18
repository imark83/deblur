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
  int nImages = 8;
  for(int k=0; k<nImages; ++k) {
    // RANDOM NUMBER GENERATOR FROM OPENCV LIB
    cv::RNG rng(1);

    // PARAMETERS
    int kernelSize = 17;
    double kernelSigma = 7;
    double sigma = 1.0e-5;
    double alpha=0.5;
    int nIter = 500;
    double mu = 1e9; //0.05 / cv::max(sigma,1.e-12);

    char fname[100];
    string imageName[] = {
      "cameraman256.png",
      "4.2.03-512.png",
      "5.3.01-1024.png",
      "brain-512.png",
      "heart-512.png",
      "lena256.png",
      "5.2.08-512.png",
      "5.3.02-1024.png"};

    cv::Mat cv_original = cv::imread("testImages/"+imageName[k], cv::IMREAD_GRAYSCALE);
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
    std::vector<double> OBJ, E, TV, S, residual;
    std::ofstream fout;

    // rop = admm1(OBJ, TV, E, S, residual, I, H, Bn, mu, alpha, nIter);
    // sprintf(fname, "testData/admm(1)-%s.txt", imageName[k].c_str());
    // std::cout << "fname = " << fname << std::endl;
    // fout.open(fname);
    // fout << "# OBJ  TV  ERROR  SNR  RESIDUAL" << std::endl;
    // for(int i=0; i<nIter; ++i)
    //   fout << i << " " << OBJ[i] << " " << TV[i]
    //         << " " << E[i] << " " << S[i] << " " << residual[i] << std::endl;
    // fout.close();


    rop = admm05(OBJ, TV, E, S, residual, I, H, Bn, mu, alpha, nIter);
    sprintf(fname, "testData/admm(05)-%s.txt", imageName[k].c_str());
    fout.open(fname);
    fout << "# OBJ  TV  ERROR  SNR  RESIDUAL" << std::endl;
    for(int i=0; i<nIter; ++i)
      fout << i << " " << OBJ[i] << " " << TV[i]
            << " " << E[i] << " " << S[i] << " " << residual[i] << std::endl;
    fout.close();


    // rop = iadmm(OBJ, TV, E, S, residual, I, H, Bn, mu, alpha, nIter);
    //
    // sprintf(fname, "testData/iadmm-%s.txt", imageName[k].c_str());
    // fout.open(fname);
    // fout << "# OBJ  TV  ERROR  SNR  RESIDUAL" << std::endl;
    // for(int i=0; i<nIter; ++i)
    //   fout << i << " " << OBJ[i] << " " << TV[i]
    //         << " " << E[i] << " " << S[i] << " " << residual[i] << std::endl;
    // fout.close();

  }

  return 0;
}
