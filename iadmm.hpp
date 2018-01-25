
#ifndef __IADMM_H__
#define __IADMM_H__

#include "matrix/cmat.hpp"

Mat iadmm(CArray &S, const Mat &img, const Mat &k, const Mat &blurred,
        double mu, double alpha, int nIter);

#endif
