#ifndef __CMAT_HPP__
#define __CMAT_HPP__

#include "mat.hpp"
#include <complex>


using namespace std;
typedef complex<double> Complex;
typedef Mat_<double> Mat;
typedef Mat_<Complex> CMat_;

class CMat : public CMat_ {
public:

  CMat();
  CMat(int r, int c) : CMat_(r, c) {}
  CMat(const CMat &op) : CMat_(op) {}
  CMat(const CMat_ &op) : CMat_(op) {}
  CMat(const Mat &op);

  CMat &operator=(const CMat_ &op);

};


// ARITHMETIC FUNCTIONS
CMat operator+(const Mat &op1, const CMat &op2);
CMat operator+(const CMat &op1, const Mat &op2);
CMat operator-(const Mat &op1, const CMat &op2);
CMat operator-(const CMat &op1, const Mat &op2);

// SCALAR PRODUCT
CMat operator*(double op1, const CMat &op2);
CMat operator*(const Complex &op1, const Mat &op2);

// MATRIX PRODUCT
CMat operator*(const CMat &op1, const Mat &op2);
CMat operator*(const Mat &op1, const CMat &op2);

// POINT PRODUCT OPERATOR
CMat operator^(const CMat &op1, const Mat &op2);
CMat operator^(const Mat &op1, const CMat &op2);

void getRand(Mat &);
void getRand(CMat &);


// COMPLEX FUNCTIONS
CMat conj(const CMat &op);
CMat abs(const CMat &op);
CMat real(const CMat &op);
double norm(const CMat &op, int n=2);
double norm2(const CMat &op);


// sign of real part
CMat sign(const CMat &op);

// max between real part of CMat and number
CMat max(const CMat &op1, double op2);

#endif
