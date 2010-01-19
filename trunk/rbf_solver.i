/* File: rbf_solver.i */
/* Interface file between python and C */

%module rbf_solver

%{
#define SWIG_FILE_WITH_INIT
extern void rbf_solver( double *x, int nx,
			double *y, int ny,
			double *g, int ng,
			double *e, int ne,
			double sigma);
%}

%include "numpy.i"
%init %{
import_array();
%}

%apply (double* IN_ARRAY1, int DIM1) {(double *x, int nx)};
%apply (double* IN_ARRAY1, int DIM1) {(double *y, int ny)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *g, int ng)};
%apply (double* IN_ARRAY1, int DIM1) {(double *e, int ne)};

extern void rbf_solver( double *x, int nx,
			double *y, int ny,
			double *g, int ng,
			double *e, int ne,
			double sigma);
%clear (double *x, int nx);
%clear (double *y, int ny);
%clear (double *g, int ng);
%clear (double *e, int ne);
