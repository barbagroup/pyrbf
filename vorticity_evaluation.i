/* File: vorticity_evaluation.i */
/* Interface file between python and C */

%module vorticity_evaluation

%{
#define SWIG_FILE_WITH_INIT
extern void vorticity_evaluation( double *wi, int nwi,
			double *xi, int nxi,
			double *yi, int nyi,
			double *xj, int nxj,
			double *yj, int nyj,
			double *gj, int ngj,
			double sigma);
%}

%include "numpy.i"
%init %{
import_array();
%}

%apply (double* INPLACE_ARRAY1, int DIM1) {(double *wi, int nwi)};
%apply (double* IN_ARRAY1, int DIM1) {(double *xi, int nxi)};
%apply (double* IN_ARRAY1, int DIM1) {(double *yi, int nyi)};
%apply (double* IN_ARRAY1, int DIM1) {(double *xj, int nxj)};
%apply (double* IN_ARRAY1, int DIM1) {(double *yj, int nyj)};
%apply (double* IN_ARRAY1, int DIM1) {(double *gj, int ngj)};

extern void vorticity_evaluation( double *wi, int nwi,
			double *xi, int nxi,
			double *yi, int nyi,
			double *xj, int nxj,
			double *yj, int nyj,
			double *gj, int ngj,
			double sigma);
%clear (double *wi, int nwi);
%clear (double *xi, int nxi);
%clear (double *yi, int nyi);
%clear (double *xj, int nxj);
%clear (double *yj, int nyj);
%clear (double *gj, int ngj);
