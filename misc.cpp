#include "misc.hpp"
#include "fft.hpp"

// FILLS WITH GAUSSIAN NOISE
void fillNoise(Mat &rop, cv::RNG &rng, double mean, double std) {
  for(int i=0; i<rop.rows; ++i) for(int j=0; j<rop.cols; ++j)
    rop(i,j) = rng.gaussian(std) + mean;

  return;
}


// CUTS VALUES OF MATRIX TO MINVAL AND MAXMAL
void mapMat(Mat &rop, double minval, double maxval) {
  for(int i=0; i<rop.rows; ++i) for(int j=0; j<rop.cols; ++j) {
    rop(i,j) =(rop(i,j) < 0.0) ? 0.0 : rop(i,j);
    rop(i,j) =(rop(i,j) > 1.0) ? 1.0 : rop(i,j);
  }
  return;
}

// GENERATES GAUSSIAN KERNEL
// size MUST be odd
Mat kernel(int size, double sigma) {
  if(size % 2 == 0) {
    std::cerr << "Kernel size must be odd" << '\n';
    exit(1);
  }
  Mat rop(size, size);
  int m = size/2;
  double d2, norm = 0.0, x, var = sigma*sigma;
  double *p =(double*) rop.data;
  for(int i=0; i<size; ++i) for(int j=0; j<size; ++j) {
    d2 =(i-m)*(i-m) +(j-m)*(j-m);
    x = exp(-0.5*d2/var);
    p[size*i+j] = x;
    norm += x;
  }
  for(int i=0; i<size; ++i) for(int j=0; j<size; ++j)
    p[size*i+j] /= norm;

  return rop;
}

// CONVOLUTES MATRIX APPLYING KERNEL USING CIRCULAR MAPPING(TORUS)
void convolute(Mat &rop, const Mat &op, const Mat &kernel) {
  // use a, b as counters for rows and cols on rop
  // use i, j as counters for rows and cols on kernel

  // COPY OP IN CASE rop == op
  Mat image(op);

  // ALLOCATE ROP AND INITIALIZE WITH ZEROS
  rop = Mat(op.rows, op.cols);
  int ii, jj;

  // MIDIUM ROW AND COL
  int mcol = kernel.cols/2;
  int mrow = kernel.rows/2;

  for(int a=0; a<rop.rows; ++a) for(int b=0; b<rop.cols; ++b) {
    // work on pixel(a,b)
    for(int i=0; i<kernel.rows; ++i) for(int j=0; j<kernel.cols; ++j) {
      // COMPUTE CIRCULAR INDECES ON OP MATRIX
      ii =(i+a-mrow);
      if(ii<0) ii+= op.rows;
      jj =(j+b-mcol);
      if(jj<0) jj+= op.cols;
      rop(a,b) += image(ii,jj)*kernel(i,j);
    }
  }
}


void psf2otf (CMat &rop, const CMat &op1, int outsize[2]) {
  // EXTEND op1 to OUTSIZE size
  int insize[2];
  insize[0] = op1.rows;
  insize[1] = op1.cols;
  rop = CMat (outsize[0], outsize[1]);

  for (int i=0; i<op1.rows; ++i) for (int j=0; j<op1.cols; ++j)
    rop(i,j) = op1(i,j);


  int toShiftR = insize[0]/2;
  int toShiftC = insize[1]/2;

  // SHIFT rows
  for (int j=0; j<rop.cols; ++j)
    for (int k=0; k<toShiftR; ++k) {
      std::complex<double> aux(rop(0,j));
      for (int i=0; i<rop.rows-1; ++i)
        rop(i,j) = rop(i+1,j);
      rop(rop.rows-1,j) = aux;
    }
  // SHIFT cols
  for (int i=0; i<rop.rows; ++i)
    for (int k=0; k<toShiftC; ++k) {
      std::complex<double> aux(rop(i,0));
      for (int j=0; j<rop.cols-1; ++j)
        rop(i,j) = rop(i,j+1);
      rop(i,rop.cols-1) = aux;
    }

  // EXECUTE FFTN
  fftn(rop);

  return;
}



// INITIALIZES CONSTANT STUFF
void getC(CMat &conjoDx, CMat &conjoDy, CMat &Nomin1, CMat &Denom1,
      CMat &Denom2, const Mat &Bn, const Mat &H) {

  int sizeB[2] = {Bn.rows, Bn.cols};

  CMat aux(1,2); aux(0,0) = 1; aux(0,1) = -1;
  CMat otfDx;
  psf2otf(otfDx, aux, sizeB);

  aux = CMat(2,1); aux(0,0) = 1; aux(1,0) = -1;
  CMat otfDy;
  psf2otf(otfDy, aux, sizeB);

  conjoDx = conj(otfDx);
  conjoDy = conj(otfDy);

  CMat otfH;
  psf2otf(otfH, H, sizeB);

  aux = CMat(Bn);
  fftn(aux);
  Nomin1 = conj(otfH)^aux;
  Denom1 = pow2(abs(otfH));
  Denom2 = pow2(abs(otfDx)) + pow2(abs(otfDy));

  return;
}



// SHRINFT FUNCTION
const double PI = 3.141592653589793;
CMat shrinft(const CMat &x, double nameta) {
  CMat rop(x.rows, x.cols);
  double aux = 3.77976314968462*pow(nameta, 0.166666666666667)/4.0;
  std::cout << "aux = " << aux << std::endl;

  std::cout << "abs(a02) = " << abs(x(0,2)) << std::endl;

  for(int i=0; i<x.rows; ++i)
    for(int j=0; j<x.cols; ++j)
      if (std::abs(x(i,j)) < aux)
        rop(i,j) = Complex(0);
      else
        rop(i,j) = (2.0/3.0)*x(i,j)*
            (1+cos(2.0*PI/3.0-2.0/3.0*
            acos(nameta/8.0*pow((abs(x(i,j))/3.0),-1.5))));
  return rop;
}


Mat toMat(const cv::Mat &op) {
  Mat rop(op.rows, op.cols);
  for(int i=0; i<rop.rows; ++i) for(int j=0; j<rop.cols; ++j) {
    rop(i,j) =((double) op.data[op.cols*i+j]) / 255.0;
  }
  return rop;
}

cv::Mat toCVMat(const Mat &op) {
  cv::Mat rop(op.rows, op.cols, CV_8UC1);
  for(int i=0; i<rop.rows; ++i) for(int j=0; j<rop.cols; ++j) {
    rop.data[rop.cols*i+j] =(int) (op(i,j) * 255.0);
  }
  return rop;
}

Mat toMat(const cv::Mat &op, double scale) {
  Mat rop(op.rows, op.cols);
  for(int i=0; i<rop.rows; ++i) for(int j=0; j<rop.cols; ++j) {
    rop(i,j) =((double) op.data[op.cols*i+j]) / scale;
  }
  return rop;
}

cv::Mat toCVMat(const Mat &op, double scale) {
  cv::Mat rop(op.rows, op.cols, CV_8UC1);
  for(int i=0; i<rop.rows; ++i) for(int j=0; j<rop.cols; ++j) {
    rop.data[rop.cols*i+j] =(int) (op(i,j) * scale);
  }
  return rop;
}


// SPECIALIZATION FOR MATLAB FUCKING VARIANCE
template <>
Complex var(Complex *op, int n) {
  Complex rop = Complex(0);
  Complex m;
  for(int i=0; i<n; ++i) {
    m = mean<Complex>(op, n);
    rop = rop + (op[i] - m) * conj(op[i]-m);
  }
  return (rop / ((double) n));
}
