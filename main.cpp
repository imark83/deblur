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



  // START WITH THE CODE
  // initialization
  double beta = 1000.0;
  CMat U(blurred);
  // int rows = blurred.rows, cols = blurred.cols;
  CMat Px;
  CMat Py;
  CMat Px0;
  CMat Py0;
  CMat PxB;
  CMat PyB;
  CMat Ux;
  CMat Uy;
  CMat Denom;
  CMat Nomin2;
  CMat Wx;
  CMat Wy;
  CMat Nono;
  CMat CImg(img);

  CMat auxX, auxY;

  CArray S(100);

  // GET CONSTANT MATRICES
  CMat conjoDx, conjoDy, Nomin1, Denom1, Denom2;
  getC (conjoDx, conjoDy, Nomin1, Denom1, Denom2, img, k);


  for(int k=0; k<200; ++k) {
    double gamma = beta / mu;
    Denom = Denom1 + gamma * Denom2;

    // initial p-problem
    Px0 = Px;
    Py0 = Py;
    // w-subproblem
    Ux = diffY(U);
    Uy = diffX(U);

    Wx = shrinft(Ux - Complex(1.0/beta)*PxB, 1.0/beta);
    Wy = shrinft(Uy - Complex(1.0/beta)*PyB, 1.0/beta);

    // u-subproblem
    auxX = Wx;
    auxY = Wy;
    fftn(auxX); fftn(auxY);
    Nomin2 = (conjoDx^auxX) + (conjoDy^auxY);

    auxX = PxB;
    auxY = PyB;
    fftn(auxX); fftn(auxY);
    Nono = (conjoDx^auxX) + (conjoDy^auxY);

    U = (Nomin1 + gamma*Nomin2 + Complex(1.0/mu)*Nono) / Denom;
    ifft(U);
    U = real(U);

    // UPDATE P
    Ux = diffY(U);
    Uy = diffX(U);

    Px = PxB + beta*(Wx-Ux);
    Py = PyB + beta*(Wy-Uy);

    PxB = Px + alpha * (Px - Px0);
    PyB = Py + alpha * (Py - Py0);

    S[k] = snr(CImg, U);

  }

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
