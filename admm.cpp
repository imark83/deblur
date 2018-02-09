#include "matrix/cmat.hpp"
#include "misc.hpp"
#include "fft.hpp"

// function [U S]= iadmm(I,H,Bn,mu,opts,alpha)

Mat admm(std::vector<double> &E, CArray &S,
        const Mat &img, const Mat &k, const Mat &blurred,
        double mu, double alpha, int nIter) {

  // START WITH THE CODE
  // initialization
  double beta = 0.01;
  CMat U(blurred);
  // int rows = blurred.rows, cols = blurred.cols;
  CMat Px(U.rows, U.cols);
  CMat Py(U.rows, U.cols);
  CMat Ux;
  CMat Uy;
  CMat Denom;
  CMat Nomin2;
  CMat Wx(U.rows, U.cols);
  CMat Wy(U.rows, U.cols);
  CMat Nono;
  CMat CImg(img);

  CMat auxX, auxY;

  S = CArray(nIter);
  E = std::vector<double>(nIter);

  // GET CONSTANT MATRICES
  CMat conjoDx, conjoDy, Nomin1, Denom1, Denom2;
  getC (conjoDx, conjoDy, Nomin1, Denom1, Denom2, blurred, k);


  for(int k=0; k<nIter; ++k) {
    std::cout << "iteracion " << k << " / " << nIter << std::endl;
    double gamma = beta / mu;
    Denom = Denom1 + gamma*Denom2;


    // w-subproblem
    Ux = diffY(U);
    Uy = diffX(U);

    Wx = sign(Ux)^max(
        abs(Ux-Px*Complex(1.0/beta))-Complex(1.0/beta),
        0.0);
    Wy = sign(Uy)^max(
        abs(Uy-Py*Complex(1.0/beta))-Complex(1.0/beta),
        0.0);

    // u-subproblem
    auxX = Wx;
    auxY = Wy;
    fftn(auxX); fftn(auxY);
    Nomin2 = (conjoDx^auxX) + (conjoDy^auxY);

    auxX = Px;
    auxY = Py;
    fftn(auxX); fftn(auxY);
    Nono = (conjoDx^auxX) + (conjoDy^auxY);

    U = (Nomin1 + gamma*Nomin2 + Complex(1.0/mu)*Nono) / Denom;
    ifftn(U);
    U = real(U);


    // UPDATE P
    Ux = diffY(U);
    Uy = diffX(U);


    Px = Px + beta*(Wx-Ux);
    Py = Py + beta*(Wy-Uy);

    S[k] = snr(CImg, U);
    E[k] = norm(CImg - U);
  }

  Mat rop(U.rows, U.cols);
  for(int i=0; i<rop.rows; ++i) for(int j=0; j<rop.cols; ++j)
    rop(i,j) = (double) std::real(U(i,j));

  return rop;

}
