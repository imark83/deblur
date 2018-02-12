#include "matrix/cmat.hpp"
#include "misc.hpp"
#include "fft.hpp"

// function [U S]= iadmm(I,H,Bn,mu,opts,alpha)

Mat iadmm(std::vector<double> &OBJ, std::vector<double> &TV,
        std::vector<double> &E, std::vector<double> &S,
        const Mat &img, const Mat &ker, const Mat &blurred,
        double mu, double alpha, int nIter) {

  // START WITH THE CODE
  // initialization
  double beta = 1000.0;
  CMat U(blurred);
  // int rows = blurred.rows, cols = blurred.cols;
  CMat Px(U.rows, U.cols);
  CMat Py(U.rows, U.cols);
  CMat Px0(U.rows, U.cols);
  CMat Py0(U.rows, U.cols);
  CMat PxB(U.rows, U.cols);
  CMat PyB(U.rows, U.cols);
  CMat Ux;
  CMat Uy;
  CMat Denom;
  CMat Nomin2;
  CMat Wx(U.rows, U.cols);
  CMat Wy(U.rows, U.cols);
  CMat Nono;
  CMat CImg(img);

  CMat auxX, auxY;
  Mat aux(U.rows, U.cols);



  S = std::vector<double>(nIter);
  E = std::vector<double>(nIter);
  TV = std::vector<double>(nIter);
  OBJ = std::vector<double>(nIter);

  // GET CONSTANT MATRICES
  CMat conjoDx, conjoDy, Nomin1, Denom1, Denom2;
  getC (conjoDx, conjoDy, Nomin1, Denom1, Denom2, blurred, ker);


  // INITIALIZE Ux, Uy
  Ux = diffY(U);
  Uy = diffX(U);

  for(int k=0; k<nIter; ++k) {
    std::cout << "iteracion " << k << " / " << nIter << std::endl;
    double gamma = beta / mu;
    Denom = Denom1 + gamma * Denom2;

    // initial p-problem
    Px0 = Px;
    Py0 = Py;
    // w-subproblem

    // Ux = diffY(U);
    // Uy = diffX(U);

    Wx = shrinft(Ux - PxB/beta, 1.0/beta);
    Wy = shrinft(Uy - PyB/beta, 1.0/beta);


    // u-subproblem
    auxX = Wx;
    auxY = Wy;
    fftn(auxX); fftn(auxY);
    Nomin2 = (conjoDx^auxX) + (conjoDy^auxY);

    auxX = PxB;
    auxY = PyB;
    fftn(auxX); fftn(auxY);
    Nono = (conjoDx^auxX) + (conjoDy^auxY);

    U = (Nomin1 + gamma*Nomin2 + Nono/mu) / Denom;
    ifftn(U);
    U = real(U);

    // UPDATE P
    Ux = diffY(U);
    Uy = diffX(U);

    Px = PxB + beta*(Wx-Ux);
    Py = PyB + beta*(Wy-Uy);

    PxB = Px + alpha * (Px - Px0);
    PyB = Py + alpha * (Py - Py0);



    // PLOT
    for(int i=0; i<aux.rows; ++i) for(int j=0; j<aux.cols; ++j)
      aux(i,j) = (double) std::real(U(i,j));
    cv::imshow("test", toCVMat(aux));
    cv::waitKey(10);

    for(int i=0; i<aux.rows; ++i) for(int j=0; j<aux.cols; ++j)
      aux(i,j) = (double) std::real(U(i,j)-CImg(i,j))*20.0;
    cv::imshow("test2", toCVMat(aux));
    cv::waitKey(10);


    // COMPUTE METADATA
    S[k] = std::real(snr(CImg, U));
    E[k] = norm(CImg - U);
    TV[k] = norm(Ux, 1) + norm(Uy, 1);

    for(int i=0; i<aux.rows; ++i) for(int j=0; j<aux.cols; ++j)
      aux(i,j) = std::real(U(i,j));
    convolute(aux, aux, ker);
    auxX = CMat(aux) - blurred;

    OBJ[k] = TV[k]
        + 0.5*mu*norm2(auxX)*norm2(auxX);
    std::cout << "tv = " << TV[k] << "   ";
    std::cout << "obj = " << OBJ[k] << std::endl << std::endl;

  }

  Mat rop(U.rows, U.cols);
  for(int i=0; i<rop.rows; ++i) for(int j=0; j<rop.cols; ++j)
    rop(i,j) = (double) std::real(U(i,j));

  return rop;

}
