
#ifndef __IADMM_H__
#define __IADMM_H__

#include "matrix/cmat.hpp"
typedef std::vector<std::complex<double> > CArray;

Mat iadmm(std::vector<double> &error, CArray &S,
        const Mat &img, const Mat &k, const Mat &blurred,
        double mu, double alpha, int nIter);

#endif
