// --------------------------------------
// SPM C++
// A. J. Smith (ajs224@cam.ac.uk)
//
//  V1.0 Extends MFA
//---------------------------------------


#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include "random.h"
#include "mfa_functions.h"
#include "mfa_params.h"
#include "Particle.h"


int main(int argc, char *argv[])
{
  
  using namespace std;
  using namespace mfaAnalytic;
  using namespace ajsRandom;

  bool loadRandState=false; // Generate a new set of seeds
  bool saveRandState=true; // Write the state so we can rewread later
 
  
  double M; // mass density
  double collProb;
  //double* xDiam = NULL;   // Pointer to double, initialize to nothing.
  int N=1024;           // Number of particles

  double lambda,lambda0;
  int alpha;
  //double t, tStop=12e0; // Current time and stopping time (default 12 seconds)
  double tau; // Waiting time

  double* rate = new double[6];  //
  int* events = new int[6];  //
  int eventsTotal;
    
  double collRate, rateTotal;
  double minProperty;

  double inFactor=10e0,outFactor=inFactor;
  
  bool firstChosen, secondChosen;
  int firstParticle, secondParticle, particleIndex;

  int moment;
  string kernelName;
  double* moments = new double[noMoments];

  double fieldMassDensCurrent,fieldMassDensFuture;

  int iter=0;
  int iterMax=1000;

  double Qin=1e0;
  double repProb;
  
  bool coagOn=true;

  bool testParticleActive=true;
    
  // Set up output file streams
  ofstream outputFile;
  ofstream momentsFile;
  ofstream diamsFile;

  string outputFileName=dataDir+"data_";
  string momentsFileName=dataDir+"moments_";
  string diamsFileName=dataDir+"diameters_";
  
  //string outputFileName="data/data_";
  //string momentsFileName="data/moments_";
  //string diamsFileName="data/diameters_";

  string ext=".txt";
  
  string desc;
  stringstream out;
  //char * pEnd;

  // Output blurb
  cout << endl;
  cout << "SPM Stochastic PBE Solver - A. J. Smith (ajs224@cam.ac.uk)" << endl;
  cout << endl;
  cout << "This code solves the continuous Smoluchowski equation with in/outflow" << endl;
  cout << "and coagulation described by constant, additive and multiplicative kernels " << endl;
  cout << "(admitting analytic solutions) in addition to a range of more physically " << endl;
  cout << "realistic kernels (run with --help for additional information)." << endl;
  cout << endl;
  
  
  // Process command line arguments
  for (int i=1; i<argc; ++i)
    {
      if (strcmp(argv[i], "--help") == 0)
	{
	  //cout << "This is the help message" << endl;
	  cout << "Usage: "<< argv[0] << " <flags>" << endl << endl;
	  cout << "where <flags> is one or more of:" << endl << endl;
	  cout << "\t" << "-t" << "\t\t" << "stopping time (default 12.0s)" << endl;
	  cout << "\t" << "-n" << "\t\t" << "number of stochastic particles (default 1024)" << endl;
	  cout << "\t" << "-in" << "\t\t" << "inflow factor (default 10.0, i.e., rate is 1/10)" << endl;
	  cout << "\t" << "-out" << "\t\t" << "outflow factor (default = inflow)" << endl;
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
	  cout << "* Still to come..."  << endl;
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
	  
	  return 0;
	}
      else if (strcmp(argv[i], "-in") == 0) 
	{
	  // Read inflow factor
	  inFactor = atof(argv[++i]); // default 1/10
	  outFactor=inFactor;
	}
      else if (strcmp(argv[i], "-out") == 0)
	{
	  // Read outflow factor
	  // If omitted inflow=outflow rate
	  outFactor = atof(argv[++i]); // default 2
	}
      else if (strcmp(argv[i], "-t") == 0)
	{
	  // Read stopping time
	  //tStop = strtod(argv[++i], &pEnd); // default 12.0s
	  tStop = atof(argv[++i]); // default 12.0s
	}
      else if (strcmp(argv[i], "-n") == 0)
	{
	  // Read numer of stochastic particles
	  N = atoi(argv[++i]); // default 1024
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

    Particle testParticle;
    Particle * fieldParticleCurrent = new Particle[N];
    Particle * fieldParticleFuture = new Particle[N];


    // Declare a Mersenne Twister random number generator
    MTRand mtrand;
    mtrand=myRand(loadRandState);
    
 

    
    // Initialise PSD to a delta delta(x-1), i.e., mono-dispersed (equivalent to setting m(x,0)=delta(x-1))
    for (int i=0; i<N; i++)
      {  
	fieldParticleCurrent[i].setMass(1e0);
	fieldParticleFuture[i].setMass(1e0);
      }

    // Initialise Phi_p
    fieldMassDensCurrent=1e0;
   
    out << N;
    desc=out.str()+"_";
    
    outputFileName+=desc+kernelName+ext;
    momentsFileName+=desc+kernelName+ext;
    diamsFileName+=desc+kernelName+ext;
    
    outputFile.open(outputFileName.c_str(), ios::out);
    momentsFile.open(momentsFileName.c_str(), ios::out);
    diamsFile.open(diamsFileName.c_str(), ios::out);
    
    //cout << "N= " << N <<endl;
    //cout << "tStop= " << tStop << "s." << endl;
    cout << "Running stochastic simulation with " << N << " stochastic particles for " << iterMax << " iterations." << endl;
    kernelName[0]=toupper(kernelName[0]); // capitalise
    cout << kernelName << " kernel selected." << endl;
    cout << "Inflow rate (1/alpha):" << 1e0/inFactor << endl;
    cout << "Outflow rate (1/beta):" << 1e0/outFactor << endl;

    cout.precision(10);
    cout.width(20);
    
    cout.setf(ios::scientific);
    outputFile.precision(8);
    momentsFile.precision(8);
    diamsFile.precision(8);
    
    // Output header
    cout << "t\t\tm0\t\t\tm1\t\t\tm2\t\t\tm3" << endl;
    
    // Iterate whilst t less than the stopping time
    while (iter<iterMax)
      {
	
	fieldMassDensFuture=0e0;

	testParticle.setMass(m_in());

	cout << "Test particle Mass=" << testParticle.getMass() << endl;
	
	while(testParticleActive)
	  {
	    	
	    //Compute rates (could shove this in a subroutine)
	    //Compute collision rate
	    collRate=0.0;
	    minProperty=1.0e8;
	    
	    for (int i=0; i<N; i++)
	      {
		minProperty=min(minProperty, fieldParticleCurrent[i].getMass());
		
		for (int j=0; j<N; j++)
		  {
		    collRate+=k(fieldParticleCurrent[i].getMass(),fieldParticleCurrent[j].getMass())/fieldParticleCurrent[j].getMass();
		  }
	      }
	    
	    collRate*=fieldMassDensCurrent/N;
	    
	    if(!coagOn)
	      collRate=0e0; // Switch off coagulation for testing purposes
	    
	    rate[0]=collRate;
	    
	    // Outflow rate
	    rate[1]=0e0;
	    for (int i=0; i<N; i++)
	      {
		rate[1]+=1e0/theta(testParticle.getMass(),outFactor);
	      }
	    
	    // Inflow rate
	    //rate[2]=N/2e0; // m_in~U[0,1] // this should always be the rate!
	    //rate[2]=N/inFactor; // m_in~U[0,1] // this should always be the rate!
	    //rate[2]=1e0/2e0;  // m_in_delta_{i1}
	    
	    // Total rate
	    //rateTotal=rate[0]+rate[1]+rate[2];
	    rateTotal=rate[0]+rate[1];
	    
	    // Compute exponentially distributed waiting time
	    tau=-log(mtrand())/(rateTotal);
	    
	    fieldMassDensFuture=fieldMassDensFuture+tau*Qin;
	    
	    repProb=eps2*Qin*tau/fieldMassDensCurrent;
	    
	    // with prob repProb replace a uniformly chosen field particle with a test particle
	    if(mtrand()<repProb)
	      {
		firstParticle=mtrand.randInt( N ) ;
		fieldParticleFuture[firstParticle].setMass(testParticle);
	      }
	    
	    // Implement binary tree, but for now, do it like a spastic
	    if(mtrand()<rate[0]/rateTotal)
	      {
		// Coagulation occurs
		// Choose a collision pair
		
		firstChosen=false;
		secondChosen=false;
	    
		while(firstChosen==false)
		  {
		    while(secondChosen==false)
		      {
			firstParticle=mtrand.randInt( N ) ;
			secondParticle=mtrand.randInt( N );
			
			//collProb=k(xDiam[firstParticle],xDiam[secondParticle])/(rateTotal*xDiam[secondParticle]);
			collProb=k(mfaParticle[firstParticle].getMass(),mfaParticle[secondParticle].getMass())/(rateTotal*mfaParticle[secondParticle].getMass());
			if( (mtrand()< collProb) && (firstParticle!=secondParticle))
			  {
			    firstChosen=true;
			    secondChosen=true; 
			  }
		      }
		    
		    
	      
		  }
		


	  mfaParticle[firstParticle]+=mfaParticle[secondParticle]; // use overloaded += operator

/*
	  cout <<"New mass is  "<< mfaParticle[firstParticle].getMass() << endl;

	  cout << "Press Enter to continue..." << endl;
	  while (1)
	  {
	      if ('\n' == getchar())
	         break;
	  }
*/
	  events[0]++;
	  testParticleActive=true;
	  
	      }
	    
	    else
	      {
		testParticleActive=false;
	      }
	  }

	// Test particle leaves the system
	
	  events[1]++;

	  //cout << "Outflow occurs" << endl;
	  
      
	  eventsTotal=events[0]+events[1]+events[2];
	  
	  cout.precision(8);
	  outputFile.precision(8);
	  momentsFile.precision(8);
	  
	  
	  
	  // Let's compute the moments of the distribution
	  for(moment=0;moment<noMoments;moment++)
	    {
	      
	      moments[moment]=0e0;
	      
	      for(int i=0;i<N;i++)
		{
		  moments[moment]+=pow(mfaParticle[i].getMass(),moment);
		}
	      moments[moment]=moments[moment]/N;
	    }
	  
      
	  //cout << tau << "\t" << t << "\t" << M << "\t" << lambda << "\t" << alpha << "\t" << events[0] << "\t" << events[1] << "\t" << events[2] << "\t" << eventsTotal << "\t" <<  firstParticle << "\t" << secondParticle << "\t" << moments[1] << endl;
	  
	  cout << t << "\t" << moments[0] << "\t" << moments[1] << "\t" << moments[2] << "\t" << moments[3]<< endl;
	  
	  
	  outputFile << tau << "\t" << t << "\t" << M << "\t" << lambda << "\t" << alpha << "\t" << events[0] << "\t" << events[1] << "\t" << events[2] << "\t" << eventsTotal << endl;
	  
	  
	  
	  momentsFile << t << "\t" << moments[0] << "\t" << moments[1] << "\t" << moments[2] << "\t" << moments[3] << "\t" << moments[4] << endl;
	  
	  
	  
      
      }

    outputFile.close();
    momentsFile.close();
    
    for(int i=0;i<N;i++)
      {
	diamsFile << mfaParticle[i].getMass() << endl;
      }
    
    diamsFile.close();
    
    //moments=computeMoments(xDiam,1);
    


    // Let's save the state of the random number generator
    if(saveRandState)
      {
	saveState(mtrand);
      }
    
    
    
    
    cout << "Events summary (collisions, inflows, outflows):" << endl << endl;
    cout << "\t" << events[0] << "\t" << events[1] << "\t" << events[2] << endl << endl;
    cout << "Simulation complete!" << endl;
    cout << "Ran stochastic simulation with " << N << " stochastic particles for " << tStop << " seconds." << endl;
    cout << kernelName << " kernel selected." << endl;
    
    // Clean up memory
    /* Leads to error.  Why?
       delete [] mfaParticle;  // When done, free memory pointed to by .
       mfaParticle = NULL;     // Clear pxDiam to prevent using invalid memory reference
    */
    delete [] rate;
    rate = NULL;
    delete [] events;
    events = NULL;
    delete [] moments;
    moments = NULL; 
    
    return 0;
    
}

