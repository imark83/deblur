#include "cmat.hpp"


CMat::CMat () {
  rows = cols = 1;
  data = new Complex[1];
  data[0] = Complex(0);
}

CMat::CMat (const Mat &op) : CMat_ (op.rows, op.cols) {
  for (int i=0; i<rows; ++i) for (int j=0; j<cols; ++j)
    (*this)(i,j) = op(i,j);
}

CMat & CMat::operator=(const CMat_ &op) {
  if(this == &op) return *this;
  if(op.rows != this->rows || op.cols != this->cols) {
    rows = op.rows; cols = op.cols;
    delete [] this->data;
    this->data = new Complex[rows*cols];
  }
  for(int i=0; i<rows; ++i) for(int j=0; j<cols; ++j)
    this->data[cols*i+j] = op.data[cols*i+j];
  return *this;
}


// COMPLEX FUNCTIONS
CMat conj(const CMat &op) {
  CMat rop(op.rows, op.cols);
  for(int i=0; i<op.rows; ++i) for(int j=0; j<op.cols; ++j)
    rop(i,j) = std::conj(op(i,j));
  return rop;
}

CMat abs(const CMat &op) {
  CMat rop(op.rows, op.cols);
  for(int i=0; i<op.rows; ++i) for(int j=0; j<op.cols; ++j)
    rop(i,j) = std::abs(op(i,j));
  return rop;
}

CMat real(const CMat &op) {
  CMat rop(op.rows, op.cols);
  for(int i=0; i<op.rows; ++i) for(int j=0; j<op.cols; ++j)
    rop(i,j) = std::real(op(i,j));
  return rop;
}

double norm(const CMat &op) {
  double rop = 0.0;
  for(int i=1; i<op.rows-1; ++i) for(int j=1; j<op.cols-1; ++j)
    rop += std::abs(op(i,j))*std::abs(op(i,j));
  return sqrt(rop);
}

// SIGN OF REAL PART
CMat sign(const CMat &op) {
  CMat rop(op.rows, op.cols);
  for(int i=0; i<op.rows; ++i) for(int j=0; j<op.cols; ++j)
    if (((double) std::real(op(i,j))) < 0) rop(i,j) = -1.0;
    else if (((double) std::real(op(i,j))) > 0) rop(i,j) = 1.0;
    else rop(i,j) = 0.0;

  return rop;
}


// max between real part of CMat and number
CMat max(const CMat &op1, double op2) {
  CMat rop(op1.rows, op1.cols);
  for(int i=0; i<rop.rows; ++i) for(int j=0; j<rop.cols; ++j)
    if((std::real(rop(i,j))) > op2)
      rop(i,j) = op1(i,j);
    else
      rop(i,j) = op2;

  return rop;
}



// OTHER FUNCTIONS
void getRand (Mat &rop) {
  for (int i=0; i<rop.rows; ++i) for (int j=0; j<rop.cols; ++j)
    rop(i,j) = (rand() % 11) - 5;
}
void getRand (CMat &rop) {
  for (int i=0; i<rop.rows; ++i) for (int j=0; j<rop.cols; ++j) {
    rop(i,j).real((rand() % 11) - 5);
    rop(i,j).imag((rand() % 11) - 5);
  }
}


// ARITHMETIC
CMat operator+ (const Mat &op1, const CMat &op2) {
  if (op1.rows != op2.rows || op1.cols != op2.cols) {
    std::cerr << "Dimensions missmatch" << '\n';
    exit(1);
  }
  CMat rop (op1.rows, op1.cols);
  for (int i=0; i<rop.rows; ++i) for (int j=0; j<rop.cols; ++j)
    rop(i,j) = op1(i,j) + op2(i,j);

  return rop;
}

CMat operator+ (const CMat &op1, const Mat &op2) {
  return op2 + op1;
}

CMat operator- (const Mat &op1, const CMat &op2) {
  if (op1.rows != op2.rows || op1.cols != op2.cols) {
    std::cerr << "Dimensions missmatch" << '\n';
    exit(1);
  }
  CMat rop (op1.rows, op1.cols);
  for (int i=0; i<rop.rows; ++i) for (int j=0; j<rop.cols; ++j)
    rop(i,j) = op1(i,j) - op2(i,j);

  return rop;
}

CMat operator- (const CMat &op1, const Mat &op2) {
  if (op1.rows != op2.rows || op1.cols != op2.cols) {
    std::cerr << "Dimensions missmatch" << '\n';
    exit(1);
  }
  CMat rop (op1.rows, op1.cols);
  for (int i=0; i<rop.rows; ++i) for (int j=0; j<rop.cols; ++j)
    rop(i,j) = op1(i,j) - op2(i,j);

  return rop;
}


// SCALAR PRODUCT
CMat operator* (double op1, const CMat &op2) {
  CMat rop (op2.rows, op2.cols);
  for (int i=0; i<op2.rows; ++i) for (int j=0; j<op2.cols; ++j)
    rop(i,j) = op1 * op2(i,j);
  return rop;
}

CMat operator* (const Complex &op1, const Mat &op2) {
  CMat rop (op2.rows, op2.cols);
  for (int i=0; i<op2.rows; ++i) for (int j=0; j<op2.cols; ++j)
    rop(i,j) = op1 * op2(i,j);
  return rop;
}

// MATRIX PRODUCT

CMat operator* (const CMat &op1, const Mat &op2) {
  if (op1.cols != op2.rows) {
    std::cerr << "Dimensions missmatch" << '\n';
    exit(1);
  }
  CMat rop (op1.rows, op2.cols);
  for (int i=0; i<rop.rows; ++i) for (int j=0; j<rop.cols; ++j) {
    rop(i,j) = 0;
    for (int k=0; k<op1.cols; ++k)
      rop(i,j) = rop(i,j) + op1(i,k)*op2(k,j);
  }
  return rop;
}

CMat operator* (const Mat &op1, const CMat &op2) {
  return op2*op1;
}

// POINT PRODUCT OPERATOR
CMat operator^ (const CMat &op1, const Mat &op2) {
  if (op1.rows != op2.rows || op1.cols != op2.cols) {
    std::cerr << "Dimensions missmatch" << '\n';
    exit(1);
  }
  CMat rop (op1.rows, op1.cols);
  for (int i=0; i<rop.rows; ++i) for (int j=0; j<rop.cols; ++j)
    rop(i,j) = op1(i,j) * op2(i,j);

  return rop;
}

CMat operator^ (const Mat &op1, const CMat &op2) {
  return op2^op1;
}
