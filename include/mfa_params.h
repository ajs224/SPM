#ifndef MFA_PARAMS_H
#define MFA_PARAMS_H

#include <vector>
#include <string>

namespace spm
{
  
  extern double Qin; // Inflow rate = \int m_in(x) dx
  extern const double eps1; // Relaxation factor 1
  extern const double eps2; // Relaxation factor 2
  
} //namespace spm


namespace mfaAnalytic
{
  enum kernelTypes { constant, additive, multiplicative, continuum, freemolecular, kinetic, shearlinear, shearnonlinear, settling, inertiasettling, berry, condensation, spmtest};

  extern kernelTypes kernelType; // Read from command line with -k flag.  Default is constant

  extern const int noMoments; // Number of moments to compute
  extern const std::string dataDir; // Output directory

} //namespace mfaAnalytic


namespace MFA
{

   // MATHEMATICAL CONSTANTS.

    extern const double PI;
    extern const double ONE_THIRD;
    extern const double TWO_THIRDS;

} //namespace MFA



#endif
