#ifndef __MAT_HPP__
#define __MAT_HPP__

#include <iostream>
#include <stdlib.h>


template <class T>
T min(T a, T b) {if(a<b) return a; return b;}
template <class T>
T max(T a, T b) {if(a>b) return a; return b;}

template <class T>
class Mat_ {
public:
  int rows;
  int cols;
  T *data;

  Mat_() {
    rows = cols = 1;
    data = new T[1];
    data[0] = T(0);
  }
  Mat_(int r, int c) {
    rows = r; cols = c;
    data = new T[r*c];
    for(int i=0; i<rows; ++i) for(int j=0; j<cols; ++j)
      data[cols*i+j] = T(0);
  }
  Mat_(const Mat_ &op) {
    rows = op.rows; cols = op.cols;
    data = new T[rows*cols];
    for(int i=0; i<rows; ++i) for(int j=0; j<cols; ++j)
      data[cols*i+j] = op.data[cols*i+j];
  }

  ~Mat_() {delete [] data;}


  // MEMBER ACCESS
  const T & operator()(int i, int j) const {
    return data[cols*i+j];
  }
  T & operator()(int i, int j) {
    return data[cols*i+j];
  }

  // ASSIGNMENT
  Mat_ & operator=(const Mat_ &op) {
    if(this == &op) return *this;
    if(op.rows != this->rows || op.cols != this->cols) {
      rows = op.rows; cols = op.cols;
      delete [] this->data;
      this->data = new T[rows*cols];
    }
    for(int i=0; i<rows; ++i) for(int j=0; j<cols; ++j)
      this->data[cols*i+j] = op.data[cols*i+j];
    return *this;
  }
};

// FRIEND OPERATORS
template <class T>
std::ostream & operator<<(std::ostream &output, const Mat_<T> &op) {
  int rows = min((int) op.rows, 6);
  int cols = min((int) op.cols, 5);
  for(int i=0; i<rows; ++i){
    for(int j=0; j<cols; ++j)
      output << op(i,j) << "  ";
    if(op.cols > 5) output << "...";
    output << std::endl;
  }
  if(op.rows > 6) output << "..." << std::endl;
  return output;
}


// ARITHMETIC OPERATORS
// ADDITION
template <class T>
Mat_<T> operator+(const Mat_<T> &op1, const Mat_<T> &op2) {
  if(op1.rows != op2.rows || op1.cols != op2.cols) {
    std::cerr << "Dimensions missmatch" << '\n';
    exit(1);
  }
  Mat_<T> rop(op1.rows, op1.cols);
  for(int i=0; i<rop.rows; ++i) for(int j=0; j<rop.cols; ++j)
    rop(i,j) = op1(i,j) + op2(i,j);

  return rop;
}

// SUBSTRACTION
template <class T>
Mat_<T> operator-(const Mat_<T> &op1, const Mat_<T> &op2) {
  if(op1.rows != op2.rows || op1.cols != op2.cols) {
    std::cerr << "Dimensions missmatch" << '\n';
    exit(1);
  }
  Mat_<T> rop(op1.rows, op1.cols);
  for(int i=0; i<rop.rows; ++i) for(int j=0; j<rop.cols; ++j)
    rop(i,j) = op1(i,j) - op2(i,j);

  return rop;
}


// SCALAR PRODUCT
template <class T>
Mat_<T> operator*(const T op1, const Mat_<T> &op2) {
  Mat_<T> rop(op2.rows, op2.cols);
  for(int i=0; i<op2.rows; ++i) for(int j=0; j<op2.cols; ++j)
    rop(i,j) = op1 * op2(i,j);
  return rop;
}
template <class T>
Mat_<T> operator*(const Mat_<T> &op1, const T op2) {
  return op2*op1;
}
// MATRIX PRODUCT
template <class T>
Mat_<T> operator*(const Mat_<T> &op1, const Mat_<T> &op2) {
  if(op1.cols != op2.rows) {
    std::cerr << "Dimensions missmatch" << '\n';
    exit(1);
  }
  Mat_<T> rop(op1.rows, op2.cols);
  for(int i=0; i<rop.rows; ++i) for(int j=0; j<rop.cols; ++j) {
    rop(i,j) = T(0);
    for(int k=0; k<op1.cols; ++k)
      rop(i,j) = rop(i,j) + op1(i,k)*op2(k,j);
  }
  return rop;
}


// POINT PRODUCT
template <class T>
Mat_<T> operator^(const Mat_<T> &op1, const Mat_<T> &op2) {
  if(op1.rows != op2.rows || op1.cols != op2.cols) {
    std::cerr << "Dimensions missmatch" << '\n';
    exit(1);
  }
  Mat_<T> rop(op1.rows, op1.cols);
  for(int i=0; i<rop.rows; ++i) for(int j=0; j<rop.cols; ++j)
    rop(i,j) = op1(i,j) * op2(i,j);

  return rop;
}

template <class T>
Mat_<T> pow2(const Mat_<T> &op) {
  Mat_<T> rop(op.rows, op.cols);
  for(int i=0; i<op.rows; ++i) for(int j=0; j<op.cols; ++j)
    rop(i,j) = op(i,j) * op(i,j);
  return rop;
}

// DIFF operator by rows. Adds final row with zeros
template <class T>
Mat_<T> diffX(const Mat_<T> &op) {
  if(op.rows < 2) {
    std::cerr << "Minimum number of rows to differenciate is 2" << '\n';
    exit(1);
  }
  Mat_<T> rop(op.rows, op.cols);
  for(int i=op.rows-2; i>=0; --i) for(int j=0; j<op.cols; ++j)
    rop(i,j) = op(i+1,j) - op(i,j);

  return rop;
}

// DIFF operator by cols. Adds final col with zeros
template <class T>
Mat_<T> diffY(const Mat_<T> &op) {
  if(op.cols < 2) {
    std::cerr << "Minimum number of cols to differenciate is 2" << '\n';
    exit(1);
  }
  Mat_<T> rop(op.rows, op.cols);
  for(int j=op.cols-2; j>=0; --j) for(int i=0; i<op.rows; ++i)
    rop(i,j) = op(i,j+1) - op(i,j);

  return rop;
}



#endif
