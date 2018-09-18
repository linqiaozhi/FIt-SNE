#ifndef NBODYFFT_H
#define NBODYFFT_H

#ifdef _WIN32
#include "winlibs/fftw3.h"
#else
#include <fftw3.h>
#endif
#include <complex>

using namespace std;

typedef double (*kernel_type)(double, double);

typedef double (*kernel_type_2d)(double, double, double, double);
	typedef double (*kerneltype)(double, double, double, double);
	typedef double (*kerneltype2d)(double, double,double,double, double, double);

    long int diff(timespec start, timespec end);
	int precompute2(double xmax, double xmin, double ymax, double ymin, int nlat, int nterms, kerneltype2d ker,double * band,double *boxl, double *boxr,  double * prods, double * xpts, double * xptsall,double *yptsall,int *irearr, fftw_complex * zkvalf );
	int nbodyfft2(int n, int ndim, double* xs, double * ys, double * charges, int nlat, int nterms,double *boxl, double *boxr,  double * prods, double * xpts, double * xptsall, double *yptsall,int* irearr, fftw_complex * zkvalf, double * pot, unsigned int nthreads);

void precompute(double y_min, double y_max, int n_boxes, int n_interpolation_points, kernel_type kernel,
                double *box_lower_bounds, double *box_upper_bounds, double *y_tilde_spacing, double *y_tilde,
                complex<double> *fft_kernel_vector);

void nbodyfft(int N, int n_terms, double *Y, double *chargesQij, int n_boxes, int n_interpolation_points,
              double *box_lower_bounds, double *box_upper_bounds, double *y_tilde_spacings, double *y_tilde,
              complex<double> *fft_kernel_vector, double *potentialsQij);

void interpolate(int n_interpolation_points, int N, const double *y_in_box, const double *y_tilde_spacings,
                 double *interpolated_values);

#endif
