
#ifndef __ADMM1_H__
#define __ADMM1_H__

#include "lazyMat/cmat.hpp"
typedef std::vector<std::complex<double> > CArray;

Mat admm1(std::vector<double> &OBJ, std::vector<double> &TV,
        std::vector<double> &E, std::vector<double> &S,
        std::vector<double> residual, const Mat &img,
        const Mat &ker, const Mat &blurred,
        double mu, double alpha, int nIter);

#endif
