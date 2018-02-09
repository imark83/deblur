
#ifndef __ADMM_H__
#define __ADMM_H__

#include "matrix/cmat.hpp"
typedef std::vector<std::complex<double> > CArray;

Mat admm(std::vector<double> &error, CArray &S,
        const Mat &img, const Mat &k, const Mat &blurred,
        double mu, double alpha, int nIter);

#endif
