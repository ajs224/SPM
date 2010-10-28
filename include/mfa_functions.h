#ifndef MFA_FUNCS_H
#define MFA_FUNCS_H

// Function headers
double k(double x, double y); // Collision kernel
double kt(double x, double y); // Majorant kernel
double* computeMoments(double* diamArray, int noMoments);
double theta(double x, double rate); // Outflow rate
double mIn(int i=1); // Inflowing particle distribution

#endif
