
#ifndef __ADMM_H__
#define __ADMM_H__

#include "matrix/cmat.hpp"
typedef std::vector<std::complex<double> > CArray;

Mat admm(std::vector<double> &OBJ, std::vector<double> &TV,
        std::vector<double> &E, std::vector<double> &S,
        const Mat &img, const Mat &ker, const Mat &blurred,
        double mu, double alpha, int nIter);

#endif
