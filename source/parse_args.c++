#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include "mfa_params.h" // Needed for kernelType

using namespace std;
using namespace mfaAnalytic;

int parseArgs(int argc, char *argv[], double & inFactor, double & outFactor, unsigned long int & iterMax, unsigned long int & N, unsigned int & L, std::string & kernelName, bool & coagOn)
{
 
  if ((argc == 1) || (strcmp(argv[1], "--help") == 0))
    {
      //cout << "This is the help message" << endl;
      cout << "Usage: "<< argv[0] << " <flags>" << endl << endl;
      cout << "where <flags> is one or more of:" << endl << endl;
      cout << "\t" << "-i" << "\t\t" << "maximum number of iterations (default 10^4)" << endl;
      cout << "\t" << "-N" << "\t\t" << "number of stochastic particles (default 1024)" << endl;
      cout << "\t" << "-L" << "\t\t" << "number of runs (default 10)" << endl;
      cout << "\t" << "-alpha" << "\t\t" << "inflow factor (default 10.0, i.e., rate is 1/10)" << endl;
      cout << "\t" << "-beta" << "\t\t" << "outflow factor (default = inflow)" << endl;
      cout << "\t" << "-k <type>" << "\t" << "kernel type, where <type> is one of:" << endl << endl;
      cout << "\t\t" << "constant" << "\t\t" << "constant kernel (default)" << endl;
      cout << "\t\t" << "additive" << "\t\t" << "additive" << endl;
      cout << "\t\t" << "multiplicative" << "\t\t" << "multiplicative" << endl;
      cout << "\t\t" << "continuum" << "\t\t" << "Brownian motion (continuum regime)" << endl;
      cout << "\t\t" << "freemolecular" << "\t\t" << "Brownian motion (free molecular regime)" << endl;
      cout << "\t\t" << "kinetic" << "\t\t\t" << "Based on kinetic theory" << endl;
      cout << "\t\t" << "shearlinear" << "\t\t" << "Shear (linear velocity profile)" << endl;
      cout << "\t\t" << "shearnonlinear" << "\t\t" << "Shear (nonlinear velocity profile)" << endl;
      cout << "\t\t" << "settling" << "\t\t" << "Gravitational settling" << endl;
      cout << "\t\t" << "inertiasettling" << "\t\t" << "Inertia and gravitational settling" << endl;
      cout << "\t\t" << "berry" << "\t\t\t" << "Analytic approximation of Berry's kernel" << endl;
      cout << "\t\t" << "condensation" << "\t\t" << "Condensation and/or branched-chain polymerisation" << endl;
      cout << "\t\t" << "spmtest" << "\t\t\t" << "Kernel used to test the Single Particle Method (SPM)" << endl;
      cout << endl;
      cout << "Examples:" << endl << endl;
      cout << "* Carry out 10 runs for 100 iterations with 512 stochastic particles, with alpha=0.1=beta and a constant kernel use:"  << endl;
      cout << "(time "<< argv[0] << " -i 100 -N 512 -L 10 -alpha 0.1 -k constant) |& tee const_coag_run_N512_L10.log" << endl << endl;
      cout << "* To do a single run for 1000 itertions with 4096 stochastic particles, with alpha=0.1=beta and the continuum kernel use:"  << endl;
      cout << "(time "<< argv[0] << " -i 1000 -N 4096 -L 1 -alpha 0.1 -k continuum) |& tee continuum_coag_run_N4096_L1.log" << endl << endl;
      
      
      /*
	cout << "* To run with inflow rate=outflow rate=1, const kernel, 16 outer loops and a max cluster size of 2^10 use:" << endl;
	cout << "\ttime "<< argv[0] << " -alpha 1 -loops 16 -p 10" << endl;
	cout << "* To run with the additive kernel and default in/outflow rates use:" << endl;
	cout << "\ttime "<< argv[0] << " -k additive -loops 256 -p 20" << endl;
	cout << "* In order to achieve convergence with more complicated kernels lower the inflow rate:" << endl;
	cout << "\ttime "<< argv[0] << " -alpha 0.1 -k multiplicative -loops 64 -p 16" << endl;
	cout << "\ttime "<< argv[0] << " -alpha 0.05 -k freemolecular -loops 64 -p 16" << endl;
      */
      
      cout << endl;
      
      return 1;
    }
  
  
  // Process command line arguments
  for (int i=1; i<argc; ++i)
    {
      if (strcmp(argv[i], "-alpha") == 0) 
	{
	  // Read inflow factor
	  inFactor = atof(argv[++i]); // default 1/10
	  outFactor=inFactor;
	}
      else if (strcmp(argv[i], "-beta") == 0)
	{
	  // Read outflow factor
	  // If omitted inflow=outflow rate
	  outFactor = atof(argv[++i]); // default 2
	}
      else if (strcmp(argv[i], "-i") == 0)
	{
	  // Read maximum number of iterations
	  //tStop = strtod(argv[++i], &pEnd); // default 12.0s
	  iterMax = atoi(argv[++i]); // default 12.0s
	}
      else if (strcmp(argv[i], "-N") == 0)
	{
	  // Read numer of stochastic particles
	  N = atoi(argv[++i]); // default 1024
	}
      else if (strcmp(argv[i], "-L") == 0)
	{
	  // Read numer of runs
	  L = atoi(argv[++i]); // default 10
	}
      else if (strcmp(argv[i], "-m") == 0)
	{
	  // print moments
	  // does this anyway
	}
      else if (strcmp(argv[i], "-k") == 0)
	{
	  // read constants appearing in multiplicative kernel k(x,y)=c*x^a*y^b
	  // c=argv[++i];
	  // a=argv[++i];
	  // b=argv[++i];

	  // Just read one of the 3 basic kernel types with analytic solution for now
	  char *kArg=argv[++i];
	  if (strcmp(kArg, "additive") == 0)
	    kernelType=additive;
	  else if (strcmp(kArg, "multiplicative") == 0)
	    kernelType=multiplicative;
	  else if (strcmp(kArg, "continuum") == 0)
	    kernelType=continuum;
	  else if (strcmp(kArg, "freemolecular") == 0)
	    kernelType=freemolecular;
	  else if (strcmp(kArg, "kinetic") == 0)
	    kernelType=kinetic;
	  else if (strcmp(kArg, "shearlinear") == 0)
	    kernelType=shearlinear;
	  else if (strcmp(kArg, "shearnonlinear") == 0)
	    kernelType=shearnonlinear;
	  else if (strcmp(kArg, "settling") == 0)
	    kernelType=settling;
	  else if (strcmp(kArg, "inertiasettling") == 0)
	    kernelType=inertiasettling;
	  else if (strcmp(kArg, "berry") == 0)
	    kernelType=berry;
	  else if (strcmp(kArg, "condensation") == 0)
	    kernelType=condensation;
	  else if (strcmp(kArg, "spmtest") == 0)
	    kernelType=spmtest;
	  else
	    kernelType=constant; // actually this is default anyway
	}
      else if (strcmp(argv[i], "-nocoag") == 0)
      	{
	  // In/outflow only
      	  coagOn=false;
      	}
    }
  
  // Find out which kernel type is selected
  switch (kernelType)
    {
    case continuum:
      kernelName="continuum";
      break;
    case freemolecular:
      kernelName="freemolecular";
      break;
    case kinetic:
      kernelName="kinetic";
      break;
    case shearlinear:
      kernelName="shearlinear";
      break;
    case shearnonlinear:
      kernelName="shearnonlinear";
      break;
    case settling:
      kernelName="settling";
      break;
    case inertiasettling:
      kernelName="inertiasettling";
      break;
    case berry:
      kernelName="berry";
      break;
    case condensation:
      kernelName="condensation";
      break;
      // Analytic kernels
    case additive:
      kernelName="additive";
      break;
    case multiplicative:
      kernelName="multiplicative";
      break;
    case spmtest:
      kernelName="spmtest";
      break;
    default: // constant kernel
      kernelName="constant";
    }

  return 0;

}
