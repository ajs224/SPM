#include "mfa_params.h"

namespace spm
{

  double Qin=1e0;
  const double eps1=2.0e-2;
  const double eps2=0.5e0;
  unsigned long int iterMax=100;
  bool coagOn=true; 
  
} //namespace spm


namespace mfaAnalytic
{
  //enum kernelTypes { constant, additive, multiplicative, continuum, freemolecular, kinetic, shearlinear, shearnonlinear, settling, inertiasettling, berry, condensation, spmtest};

  kernelTypes kernelType; // Read from command line with -k flag.  Default is constant

  const int noMoments=4; // Number of moments to compute
  const std::string dataDir="data/"; // Output directory

} //namespace mfaAnalytic


namespace MFA
{

   // MATHEMATICAL CONSTANTS.

    const double PI = 3.1415926535897932384626433832795;
    const double ONE_THIRD  = 3.3333334e-01;
    const double TWO_THIRDS = 6.6666667e-01;

} //namespace MFA
