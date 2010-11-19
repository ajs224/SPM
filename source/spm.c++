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
#include <map>
#include <vector>
#include "random.h"
#include "mfa_functions.h"
#include "mfa_params.h"
#include "Particle.h"


int main(int argc, char *argv[])
{
  
  using namespace std;
  using namespace spm;
  using namespace mfaAnalytic;
  using namespace ajsRandom;

  bool loadRandState=false; // Generate a new set of seeds
  bool saveRandState=true; // Write the state so we can rewread later
   
  //double M; // mass density
  double collProb;
  //double* xDiam = NULL; // Pointer to double, initialize to nothing.
  unsigned long int N=1024; // Number of field particles

  //double lambda,lambda0;
  //int alpha;
  double t, tau; // Current time and waiting time

  double* rate = new double[6]; //
  //int* events = new int[6]; //
  map<string,unsigned long int> event;

  vector<string> eventTypes;
    
  eventTypes.push_back("coagulation");
  eventTypes.push_back("inflow");
  eventTypes.push_back("outflow");

  //string kernelLabels[]={"constant", "additive"};

  unsigned long int eventsTotal = 0;
    
  double collRate, rateTotal;
  double minProperty;

  double inFactor=2e0,outFactor=inFactor;
  
  bool particleChosen;
  unsigned long int particleIndex;
  
  int moment;
  string kernelName;
  double* moments = new double[noMoments];
  
  double fieldMassDensCurrent,fieldMassDensFuture;
  
  unsigned long int iter=1;
  //int iterMax=1000;
  
  //double Qin=1e0;
  double repProb;
  
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
	  cout << "\t" << "-i" << "\t\t" << "maximum number of iterations (default 10^4)" << endl;
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
      else if (strcmp(argv[i], "-i") == 0)
	{
	  // Read iterMax
	  //tStop = strtod(argv[++i], &pEnd); // default 12.0s
	  iterMax = atoi(argv[++i]); // default 12.0s
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

    // Declare particle types used in the stochastic simulation
    Particle testParticle;
    Particle * fieldParticleCurrent = new Particle[N];
    Particle * fieldParticleFuture = new Particle[N];
    
    // Declare a Mersenne Twister random number generator
    MTRand mtrand;
    mtrand=myRand(loadRandState);
    
    // Initialise PSD to a delta delta(x-1), i.e., mono-dispersed (equivalent to setting m(x,0)=delta(x-1))
    for (unsigned long int i=0; i<N; i++)
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
    cout << "Running SPM stochastic simulation with " << N << " stochastic particles for " << iterMax << " iterations." << endl;
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

    
    // Populate event map event["<event type>"]=# events of this type
    for(vector<string>::iterator it=eventTypes.begin();it!=eventTypes.end();++it)
      event[*it] = 0;
    
    // Iterate whilst t less than the stopping time
    while (iter<iterMax)
      {
	
	// Initialise time
	t=0e0;
	
	fieldMassDensFuture=0e0;

	// Generate a test particle sampled from m_in
	//testParticle.setMass(mIn());
	testParticle.setMass(mtrand());
	//cout << "Test particle Mass=" << testParticle.getMass() << endl;

	// Update events map
	event["inflow"]++;
	
	while(testParticleActive)
	  {
	    	
	    //Compute rates (could shove this in a subroutine)
	    //Compute collision rate
	    collRate=0.0;
	    minProperty=1.0e8;
	    
	    for (unsigned long int i=0; i<N; i++)
	      {
		minProperty=min(minProperty, fieldParticleCurrent[i].getMass());
		collRate+=k(testParticle.getMass(),fieldParticleCurrent[i].getMass())/fieldParticleCurrent[i].getMass();
		
	      }
	    
	    collRate*=fieldMassDensCurrent/N;
	    
	    if(!coagOn)
	      collRate=0e0; // Switch off coagulation for testing purposes
	    
	    rate[0]=collRate;
	    
	    // Outflow rate
	    
	    // To match f90 SPM (update this once fully working)
	    rate[1]=0.1*testParticle.getMass();
	    
	    
	    /*
	      rate[1]=0e0;
	      for (unsigned long int i=0; i<N; i++)
	      {
	      rate[1]+=1e0/theta(testParticle.getMass(),outFactor);
	      }
	      /*
	      
	      
	    // Inflow rate
	    //rate[2]=N/2e0; // m_in~U[0,1] // this should always be the rate!
	    //rate[2]=N/inFactor; // m_in~U[0,1] // this should always be the rate!
	    //rate[2]=1e0/2e0;  // m_in_delta_{i1}
	    
	    // Total rate
	    //rateTotal=rate[0]+rate[1]+rate[2];
	    rateTotal=rate[0]+rate[1];
	    
	    // Compute exponentially distributed waiting time
	    tau=-log(mtrand())/(rateTotal);
	    
	    // Wait time tau
	    t+=tau;
	    
	    // Recalculate future field particles mass density
	    fieldMassDensFuture+=tau*Qin;

	    /*
	    cout << "fieldMassDensCurrent=" << fieldMassDensCurrent << " fieldMassDensFuture=" << fieldMassDensFuture << endl;
	    cout << "tau=" << tau << " Qin=" << Qin << " rateTotal=" << rateTotal << endl;
	    
	    cin.sync();
	    cin.get();
	    */
	    
	    repProb=eps2*Qin*tau/fieldMassDensCurrent;

	    // Check the rates
	    //cout << rate[0] << " " << rate[1] << " " << rateTotal << " " << tau << " " << t << " " << repProb << " " << Qin << " " << rate[0]/rateTotal << endl;
	    

	    // With probability repProb replace a uniformly chosen field particle with a test particle
	    if(mtrand()<repProb)
	      {
		particleIndex=mtrand.randInt( N ) ;
		fieldParticleFuture[particleIndex].setMass(testParticle.getMass());
		//fieldParticleFuture[firstParticle]=testParticle;
	      }

	    // Particle collision
	    // The test particle is active whilst it continues to collide
	    // it become inactive once no more collision occur and it leaves the domain

	    /*
	      double tst=mtrand();
	      double pp=rate[0]/rateTotal;
	      
	      cout << "tst=" << tst << " pp=" << pp << endl << endl;
	      cin.sync();
	      cin.get();
	    */
	    
	    // Implement binary tree, but for now, do it like a spastic
	    if(mtrand()<rate[0]/rateTotal)
	      {

		//cout << "Coagulation occurs!" << endl << endl;
		
		
		// Coagulation occurs
		// Choose a collision pair
		
		particleChosen=false;
		
		while(particleChosen==false)
		  {
		    
		    particleIndex=mtrand.randInt( N );
		    
		    //collProb=k(testParticle.getMass(), fieldParticleCurrent[particleIndex].getMass())* \
		    //  fieldMassDensCurrent/(N*rateTotal*fieldParticleCurrent[particleIndex].getMass());

		    double kMaj=k(testParticle.getMass(), minProperty)/minProperty;

		    collProb=(k(testParticle.getMass(), fieldParticleCurrent[particleIndex].getMass())/fieldParticleCurrent[particleIndex].getMass())/kMaj;
		    double tmpRand=mtrand();
		    
		    
		    if(tmpRand<collProb)
		      particleChosen=true;

		    
		    cout << "tp mass = " << testParticle.getMass() <<  " particleIndex = " << particleIndex << " tmpRand=" << tmpRand << endl;
		    cout << "collProb = " << collProb << " rateTotal=" << rateTotal << " N=" << N << endl;
		    cout << "k=" << k(testParticle.getMass(), fieldParticleCurrent[particleIndex].getMass()) << endl;
		    cout << "fieldMassDensCurrent=" << fieldMassDensCurrent << " y_{p,i}=" << fieldParticleCurrent[particleIndex].getMass() << endl;
		    cout << "minProperty=" << minProperty << endl;
		    cout << "kMaj=" << kMaj << endl;
		    cin.sync();
		    cin.get();
		  }

		
		//cout <<"#coags= "<< event["coagulation"] << " tp mass = " << testParticle.getMass() <<  " collProb = " << collProb << endl;
		//testParticle.setMass(testParticle.getMass()+fieldParticleCurrent[particleIndex].getMass());
		testParticle+=fieldParticleCurrent[particleIndex]; // use overloaded += operator

		testParticleActive=true; // collision occurred, so test particle still active
		
		// Update events map
		event["coagulation"]++;
	      }
	    else
	      {
		testParticleActive=false;
	      }

	    //cout << "tp mass = " << testParticle.getMass() << " collProb = " << collProb << "rate[0]/rateTotal = " << rate[0] << endl;
	  
	    
	    //cout << "TP active:" << testParticleActive << " tp mass = " << testParticle.getMass() << " collProb = " << collProb << endl;
   
	  }
	
	//cout << "fieldMassDensCurrent=" << fieldMassDensCurrent << " fieldMassDensFuture=" << fieldMassDensFuture <<endl;
	
	// Test particle leaves the system

	// Underrelax the field particles mass density
	fieldMassDensCurrent=eps1*fieldMassDensFuture+(1e0-eps1)*fieldMassDensCurrent;
	
	//cout << "fieldMassDensCurrent=" << fieldMassDensCurrent << " fieldMassDensFuture=" << fieldMassDensFuture <<endl;
	//cin.sync();
	//cin.get();
	

	
	
	// Update events map
	event["outflow"]++;

	/*
	for(vector<string>::iterator it=eventTypes.begin();it!=eventTypes.end();++it)
	  cout << *it << " ";
	for(vector<string>::iterator it=eventTypes.begin();it!=eventTypes.end();++it)
	  cout << "\t " << event[*it];
	cout << endl;
	*/

	
	//cout << "Outflow occurs" << endl;
	
	// Update current field particles array
	for(unsigned long int i=0; i<N; i++)
	  {  
	    fieldParticleCurrent[i]=fieldParticleFuture[i]; // overloaded = op
	  }
	
    	// Add up the total number of events of each type
	for(map<string,unsigned long int>::const_iterator it=event.begin();it!=event.end();++it)
	  eventsTotal+=it->second;

	
	cout.precision(8);
	outputFile.precision(8);
	momentsFile.precision(8);
	
	// Let's compute the moments of the distribution
	for(moment=0;moment<noMoments;moment++)
	  {
	    
	    moments[moment]=0e0;
	    
	    for(unsigned long int i=0;i<N;i++)
	      {
		moments[moment]+=pow(fieldParticleCurrent[i].getMass(),moment);
	      }
	    moments[moment]=moments[moment]/N;
	  }
	
	
	//cout << iter << "\t" << t << "\t" << moments[0] << "\t" << moments[1] << "\t" << moments[2] << "\t" << moments[3]<< endl;

	cout << iter << "\t" << event["coagulation"] << "\t" << fieldMassDensCurrent << endl;
		
	//outputFile << tau << "\t" << t << "\t" << M << "\t" << alpha << "\t" << events[0] << "\t" << events[1] << "\t" << events[2] << "\t" << eventsTotal << endl;
	momentsFile << t << "\t" << moments[0] << "\t" << moments[1] << "\t" << moments[2] << "\t" << moments[3] << "\t" << moments[4] << endl;

	// Update iteration counter
	iter++;
	
      } // End of iteration
    
    
    outputFile.close();
    momentsFile.close();
    
    for(unsigned long int i=0;i<N;i++)
      {
	diamsFile << fieldParticleCurrent[i].getMass() << endl;
      }
    
    diamsFile.close();
    
   
    // Let's save the state of the random number generator
    if(saveRandState)
      {
	saveState(mtrand);
      }  
    
    cout << "Events summary (";
    for(vector<string>::iterator it=eventTypes.begin();it!=eventTypes.end();++it)
      cout << *it << " ";
    cout << "):" << endl << endl;
    for(vector<string>::iterator it=eventTypes.begin();it!=eventTypes.end();++it)
      cout << "\t " << event[*it];
    cout << endl << endl;
    
    
    cout << "Simulation complete!" << endl;
    cout << "Ran stochastic simulation with " << N << " stochastic particles for " << iterMax << " iterations." << endl;
    cout << kernelName << " kernel selected." << endl;
    
    // Clean up memory
   
    delete [] rate;
    rate = NULL;
    //delete [] events;
    // events = NULL;
    delete [] moments;
    moments = NULL; 
    
    return 0;
    
}

