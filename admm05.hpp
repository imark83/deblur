
#ifndef __ADMM05_H__
#define __ADMM05_H__

#include "matrix/cmat.hpp"
typedef std::vector<std::complex<double> > CArray;

Mat admm05(std::vector<double> &OBJ, std::vector<double> &TV,
        std::vector<double> &E, std::vector<double> &S,
        std::vector<double> &residual, const Mat &img,
        const Mat &ker, const Mat &blurred,
        double mu, double alpha, int nIter);

#endif
