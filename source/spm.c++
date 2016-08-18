// ----------------------------------------------
// MFA C++
// A. J. Smith (ajs224@cam.ac.uk)
//
// Multiple runs 
//-----------------------------------------------

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>

#include "random.h"
#include "mfa_functions.h"
#include "mfa_params.h"
#include "Particle.h"

using namespace spm;
int main(int argc, char *argv[])
{
  
  using namespace std;
  using namespace spm;
  using namespace mfaAnalytic;
  using namespace ajsRandom;
  
  bool loadRandState=false; // Generate a new set of seeds
  bool saveRandState=true; // Write the state so we can rewread later
  
  // Declare a Mersenne Twister random number generator
  MTRand mtrand;
  
  double M; // mass density
  double collProb;
  //double* xDiam = NULL;   // Pointer to double, initialize to nothing.
  unsigned long int N = 1024;           // Number of particles
  unsigned int L = 10;// Default number of runs
     
  unsigned int noProcesses = 6; 
  double lambda,lambda0;
  int alpha;
  double t, tStop=12e0; // Current time and stopping time (default 12 seconds)
  double tau; // Waiting time
  
  unsigned long int eventsTotal = 0;
  
  double collRate, inRate, outRate, rateTotal;
  double minProperty;
  
  //double inFactor=2e0,outFactor=inFactor;
  double inFactor;
  double outFactor;

  
  bool particleChosen;
  bool firstChosen, secondChosen;
  int firstParticle, secondParticle, particleIndex;

  string kernelName;
  
  bool coagOn=true;


  // SPM specific - Could use MFA inheritence here
  double fieldMassDensCurrent,fieldMassDensFuture;
  unsigned long int iter=1;
  unsigned long int iterMax = 10000;
  double repProb;
  bool testParticleActive=true;
   
  // Set up output file streams
  ofstream outputFile;
  ofstream momentsFile;
  ofstream diamsFile;

  string ext=".txt";
  
  string desc;
  stringstream out;
  
  // Set output precision
  cout.precision(10);
  cout.width(20);
  cout.setf(ios::scientific);

  momentsFile.precision(10);
  momentsFile.width(20);
  momentsFile.setf(ios::scientific);
  
  outputFile.precision(8);
  diamsFile.precision(8);

  // Print program header blurb
  cout << endl;
  cout << "SPM Stochastic PBE Solver - A. J. Smith (ajs224@cam.ac.uk)" << endl;
  cout << endl;
  cout << "This code solves the continuous Smoluchowski equation with in/outflow" << endl;
  cout << "and coagulation described by constant, additive and multiplicative kernels " << endl;
  cout << "(admitting analytic solutions) in addition to a range of more physically " << endl;
  cout << "realistic kernels (run with --help for additional information)." << endl;
  cout << endl;

  // Parse command line arguments
  if(parseArgs(argc, argv, inFactor, outFactor, iterMax, N, L, kernelName, coagOn))
    {
      // Didn't enter any arguments
      return 0;
    }
  
  // Setup output file names
  string outputFileName=dataDir+kernelName+"_data";
  string momentsFileName=dataDir+kernelName+"_moments";
  string diamsFileName=dataDir+kernelName+"_diameters";
  
  out << "_alpha" << inFactor << "_beta" << outFactor << "_N" << N << "_L" << L;
  desc=out.str()+"";
  
  outputFileName += desc + ext;
  momentsFileName += desc + ext;
  diamsFileName += desc + ext;
  
  // Open file handles
  outputFile.open(outputFileName.c_str(), ios::out);
  momentsFile.open(momentsFileName.c_str(), ios::out);
  diamsFile.open(diamsFileName.c_str(), ios::out);
  
  // Declare a Mersenne Twister random number generator
  //MTRand mtrand;
  mtrand=myRand(loadRandState);
  
  // 3D array containing moments at all timesteps for all runs
  vector<vector<vector<double> > > momentsAllRuns(L, vector<vector<double> >(noMoments+1));
  vector<vector<vector<double> > > momentsAllRunsInterp(L, vector<vector<double> >(noMoments)); // don't need to store time here because they're all the same
  vector<double> interpTimes;
  
  vector<double> moments(noMoments,0e0); // vector moments at current timestep of current run  
   
  vector<double> rate(noProcesses,0e0); // vector of rates 
  vector<int> events(noProcesses,0);  // vector to keep track of number of each event/process which have occurred

  map<string,unsigned long int> event;
  
  vector<string> eventTypes;
  
  eventTypes.push_back("coagulation");
  eventTypes.push_back("inflow");
  eventTypes.push_back("outflow");

  // SPM Specific vectors
  Particle testParticle;
  vector<Particle> fieldParticleCurrent(N,Particle()); // vector of stochastic particles
  vector<Particle> fieldParticleFuture(N,Particle()); // vector of stochastic particles
  
  cout << "Running stochastic SPM simulation with " << N << " stochastic field particles for " << iterMax << " iterations." << endl;
  kernelName[0]=toupper(kernelName[0]); // capitalise
  cout << kernelName << " kernel selected." << endl;
  cout << "Inflow rate (1/alpha):" << 1e0/inFactor << endl;
  cout << "Outflow rate (1/beta):" << 1e0/outFactor << endl;
    
  // Output moments file header
  momentsFile << "t \t m0 \t m0_err\t m1 \t m1_err\t m2 \t m2_err\t m3 \t m3_err\t m4 \t m4_err" << endl;

  // Start the simulation
  for(int run = 0; run < L; run++)
    {
      
      t=0e0;
      
      lambda0=1e0/N;
      
      alpha=0;
            
      // Initialise PSD to a delta delta(x-1), i.e., mono-dispersed (equivalent to setting m(x,0)=delta(x-1))
      // fancy method
      /*
     for (std::pair<iterCurrent, iterFuture> fieldParticle(fieldParticleCurrent.begin(), fieldParticleFuture.begin()); fieldParticle.first != fieldParticleCurrent.end(); ++fieldParticle.first, ++fieldParticle.second)
	{
	  fieldParticle->first->setIndex(iterCurrent-fieldParticleCurrent.begin());
	  fieldParticle->first->set((double) 1e0);
	  
	  fieldParticle->second->setIndex(iterFuture-fieldParticleFuture.begin());
	  fieldParticle->second->set((double) 1e0);
	}
	*/

      // Initialise PSD to a delta delta(x-1), i.e., mono-dispersed (equivalent to setting m(x,0)=delta(x-1))
      // Less fancy method
      for (unsigned long int i=0; i<N; i++)
	{  
	  fieldParticleCurrent[i].setMass(1e0);
	  fieldParticleFuture[i].setMass(1e0);
	}
      
      // Initialise Phi_p
      fieldMassDensCurrent=1e0;
      
      // Populate event map event["<event type>"]=# events of this type
      for(vector<string>::iterator it=eventTypes.begin();it!=eventTypes.end();++it)
	event[*it] = 0;
      
      cout << "Performing run " << run + 1 << " of " << L << endl;
      
      // Output header
      cout << "t\t\tm0\t\t\tm1\t\t\tm2\t\t\tm3" << endl;
      




      // Iterate whilst t less than the stopping time
      while (iter<iterMax)
	{
	  
	  // Initialise time
	  t=0e0;
	  
	  fieldMassDensFuture=0e0;

	  // Generate a test particle sampled from m_in
	  //testParticle.setMass(mIn());
	  testParticle.setMass(mtrand());  // currently ~ U[0,1]
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
	      */
	      
	      
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
		      
		      if(false)
			{
			  
			  cout << "tp mass = " << testParticle.getMass() <<  " particleIndex = " << particleIndex << " tmpRand=" << tmpRand << endl;
			  cout << "collProb = " << collProb << " rateTotal=" << rateTotal << " N=" << N << endl;
			  cout << "k=" << k(testParticle.getMass(), fieldParticleCurrent[particleIndex].getMass()) << endl;
			  cout << "fieldMassDensCurrent=" << fieldMassDensCurrent << " y_{p,i}=" << fieldParticleCurrent[particleIndex].getMass() << endl;
			  cout << "minProperty=" << minProperty << endl;
			  cout << "kMaj=" << kMaj << endl;
			  cin.sync();
			  cin.get();
			}
		      
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

	  
	  // Let's compute the moments of the distribution
	  for(int moment=0;moment<noMoments;moment++)
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
	  
	} // End of single run
    



      /*
      // MFA CODE
      // Iterate whilst t less than the stopping time
      while (t<tStop)
	{
	  
	  lambda=lambda0*pow(double( N )/(N-1), alpha);
	  
	  M=lambda*N;
	  
	  // Populate momentsAllRuns array with current timestep
	  momentsAllRuns[run][0].push_back(t);
	  
	  // Let's compute the moments of the distribution
	  for(vector<double>::iterator iterMom=moments.begin();iterMom != moments.end();++iterMom)
	    {	       
	      *iterMom=0e0;
	      int moment = iterMom - moments.begin();
	      
	      for(vector<Particle>::iterator iterPart = mfaParticle.begin(); iterPart != mfaParticle.end(); ++iterPart)
		{
		  *iterMom+=pow(iterPart->getMass(), moment);
		}
	      
	      *iterMom /= N;
	      
	      // Populate momentsAllRuns array with moments for current timestep and current runs
	      momentsAllRuns[run][moment+1].push_back(*iterMom);
	      
	    }
	  
	  
	  // Output data to screen and file
	  cout << t << "\t";
	  for(vector<double>::iterator iterMom=moments.begin();iterMom != moments.end();++iterMom)
	    cout << *iterMom << "\t";
	  cout << endl;
	  
	  outputFile << run <<  "\t" << tau
		     << "\t" << t << "\t"
		     << M << "\t" << lambda << "\t"
		     << alpha << "\t" << events[0]
		     << "\t" << events[1] << "\t"
		     << events[2] << "\t" << eventsTotal << endl;
	  
	  
	  //Compute rates (could shove this in a subroutine)
	  //Compute collision rate
	  collRate=0.0;
	  minProperty=1.0e8;
	  
	  for(vector<Particle>::iterator iterPart1 = mfaParticle.begin(); iterPart1 != mfaParticle.end(); ++iterPart1)
	    {
	      minProperty=min(minProperty, iterPart1->getMass());
	      
	      for(vector<Particle>::iterator iterPart2 = mfaParticle.begin(); iterPart2 != mfaParticle.end(); ++iterPart2)
		{
		  collRate+=k(iterPart1->getMass(),iterPart2->getMass())/iterPart2->getMass();
		}
	    }
	  
	  collRate*=lambda;
	  
	  if(!coagOn)
	    collRate=0e0; // Switch off coagulation for testing purposes
	  
	  // Outflow rate
	  outRate=0e0;
	  for(vector<Particle>::iterator iterPart = mfaParticle.begin(); iterPart != mfaParticle.end(); ++iterPart)
	    {
	      outRate += 1e0/theta(iterPart->getMass(),outFactor);
	    }
	  
	  // Inflow rate
	  inRate = N / inFactor; // m_in~U[0,1] // this should always be the rate!
	  //inRate=N/2e0; // m_in~U[0,1] // this should always be the rate!
	  //inRate = 1e0/2e0;  // m_in_delta_{i1}
	  
	  rate[0] = collRate;
	  rate[1] = outRate;
	  rate[2] = inRate;
	  
	  // Total rate
	  rateTotal=rate[0]+rate[1]+rate[2];
	  
	  // Compute exponentially distributed waiting time
	  tau=-log(mtrand())/(rateTotal);
	  
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
		      firstParticle=mtrand.randInt( N-1 ) ;
		      secondParticle=mtrand.randInt( N-1 );
		      
		      collProb=k(mfaParticle[firstParticle].getMass(),mfaParticle[secondParticle].getMass())/
			(rateTotal*mfaParticle[secondParticle].getMass());
		      
		      if( (mtrand()< collProb) && (firstParticle != secondParticle))
			{
			  firstChosen=true;
			  secondChosen=true; 
			}
		    }		    
		}
	      
	      
	      mfaParticle[firstParticle] += mfaParticle[secondParticle]; // use overloaded += operator
	      
	      events[0]++;	  
	      
	    } // Coagulation process complete 
	  
	  else if( mtrand() < (rate[1]/(rate[1]+rate[2])) )
	    {
	      
	      // Particle leaves the system
	      firstChosen=false;
	      secondChosen=false;
	      
	      while(firstChosen == false)
		{
		  
		  firstParticle=mtrand.randInt( N-1 );
		  // should select according to outflow distribution
		  firstChosen=true;
		}
	      
	      while(secondChosen == false)
		{
		  secondParticle=mtrand.randInt( N-1 );
		  if(firstParticle != secondParticle)
		    secondChosen=true;
		}
	      
	      // Replace particle i with a couple of particle j (that is, particle i leaves the system)
	      //mfaParticle[firstParticle].getMass()=mfaParticle[secondParticle].getMass();
	      mfaParticle[firstParticle] = mfaParticle[secondParticle]; // This ought to do a shallow copy, copying all the member data
	      alpha--;
	      events[1]++;
	      
	      //cout << "Outflow occurs" << endl;
	      
	    } // Outflow process complete
	  else
	    {
	      // Particle inflow
	      particleIndex=mtrand.randInt( N-1 );
	      //mfaParticle[particleIndex].setMass(mtrand()); // m_in~U[0,1]
	      mfaParticle[particleIndex].setMass(mIn()); // m_in_delta_{i1}
	      alpha++;
	      events[2]++;
	      
	      //cout << "Inflow occurs" << endl;
	      
	    } // Inflow process complete
	  
	  // Update stochastic time
	  t+=tau;
	  
	  eventsTotal=events[0]+events[1]+events[2];
	  
	  
	}// End of current run	   
      */      
      
      
      
      /*
      cout << endl << "PSD:" << endl;
      for(vector<Particle>::iterator iterPart = mfaParticle.begin(); iterPart != mfaParticle.end(); ++iterPart)
	{
	  cout  << "Particle "<< iterPart- mfaParticle.begin() 
	  << " has index " <<  iterPart->getIndex()<< " has mass " << iterPart->getMass() << endl;
	  
	}
      */

      // Print summary of events for current run
      cout << "Events summary (collisions, inflows, outflows):" << endl << endl;
      cout << "\t" << events[0] << "\t" << events[1] << "\t" << events[2] << endl << endl;
      cout << "Run " << run + 1 << " of "<< L << " complete!" << endl;
      
    } // End of L runs
  

  return 0;
  

  // Post-process momentsAllRuns array here
  // In order to average over all the runs we need to interpolate the points (different runs have different #s of timesteps)
  // Strategy is to assume a linear time spacing with a number of steps equal to the # of steps in the largest run

  cout << "Post-processing data..." << endl;
  
  int maxNoSteps = 0;
  for(vector<vector<vector<double> > >::iterator iterRun = momentsAllRuns.begin(); iterRun != momentsAllRuns.end(); ++iterRun)
    {
      if ((*iterRun)[0].size() > maxNoSteps)
	maxNoSteps = (*iterRun)[0].size();
    }
  //cout << "Max number of steps is " << maxNoSteps << endl;
  
  t = 0;
  
  int lStep, rStep;
  double t1, t2;
  double interpedMom;
  
  /* - Test interpolation routine
     cout << "interpMon(1.5, 1, 2, 2, 4) = " << interpMon(1.5, 1, 2, 2, 4) << endl;
     cout << "interpMon(1, 1, 2, 2, 4) = " << interpMon(1, 1, 2, 2, 4) << endl;
     cout << "interpMon(2, 1, 2, 2, 4) = " << interpMon(2, 1, 2, 2, 4) << endl;
     cout << "interpMon(1.75, 1, 2, 2, 4) = " << interpMon(1.75, 1, 2, 2, 4) << endl;
     return 0 ;
  */

  cout << "Interpolating data... ";
  
  for (int timeStep = 0; timeStep < maxNoSteps; timeStep++)
    {
      interpTimes.push_back(t);
      t += (double) tStop/maxNoSteps;
      //cout << "interpolation t = " << t << endl;
      
      // Need to loop through all runs and find which steps t is between
      
      for (int run = 0; run < L; run++)
	{
	  lStep = 0;
	  rStep = 0;
	  for (int runStep = 0; runStep < momentsAllRuns[run][0].size() - 1; runStep++)
	    {
	      // We're trying to find the two steps in the current run in which t lies
	      if( (t>momentsAllRuns[run][0][runStep]) && (t<momentsAllRuns[run][0][runStep+1]) )
		{
		  lStep = runStep;
		  rStep = runStep+1;
		}
	      else if(t > momentsAllRuns[run][0][runStep+1])
		{
		  lStep = runStep+1;
		  rStep = runStep+1;
		}
	    }
	  //cout << "run " << run << " first t = " << momentsAllRuns[run][0][1];
	  //cout << " # timesteps = " << momentsAllRuns[run][0].size() << endl;
	  
	  if (lStep == 0 && rStep == 0)
	    {
	      cout << "Interpolation failed!" << endl;
	      cout << "lStep" << lStep << endl;
	      cout << "rStep" << rStep << endl;
	      //cout << "" << << endl;
	      
	      return 0;
	    }
	  
	  t1 = momentsAllRuns[run][0][lStep];
	  t2 = momentsAllRuns[run][0][rStep];
	  
	  
	  //momentsAllRunsInterp[run][0][timeStep + 1] = t;
	  //momentsAllRunsInterp[run][0].push_back(t);
	  
	  for (int mom = 1; mom <= noMoments; mom++)
	    {
	      if (lStep != rStep)
		interpedMom = interpMon(t,t1, t2, momentsAllRuns[run][mom][lStep],momentsAllRuns[run][mom][rStep]);
	      else
		interpedMom = momentsAllRuns[run][mom][rStep];
	      
	      //momentsAllRunsInterp[run][mom][timeStep + 1] = interpedMom;
	      momentsAllRunsInterp[run][mom-1].push_back(interpedMom); // -1 because we don't store time now
	    }
	  
	  //cout << "t = " << t << " lies between " << momentsAllRuns[run][0][lStep] << " and " << momentsAllRuns[run][0][rStep]
	  //     << "steps " << lStep << " and " << rStep << endl;
	  
	}    
    }
  
  cout << "Done!" << endl;

  /*
    for(vector<double>::iterator iterT= interpTimes.begin(); iterT != interpTimes.end(); ++iterT)
     {
     cout << " interpolated time = " << *iterT << endl;
     
     }
  */

  cout << "Averaging moment data...";
  
  // t = 0 are the ICs, so no need to interpolate here, just output values from the first run
  // and zero variances (since ICs are identical for all runs)
  momentsFile << "0.0000000000e+00\t";
  for(int mom = 0; mom < noMoments; mom++)
    momentsFile << momentsAllRuns[0][mom+1][0] << "\t0.0000000000e+00\t";
  
  momentsFile << endl;  

  double momentCum;
  double momentSqrdCum;
  double momMean, momVar;
  
  
  // We now have an array momentsAllRunsInterp of interpolated moments at the *same* time, so we can loop
  // through and calculate confidence intervals
  
  for (int timeStep = 0; timeStep < maxNoSteps; timeStep++)
    {
      momentsFile << interpTimes[timeStep] << "\t";
      for (int mom = 0; mom < noMoments; mom++)
	{
	  momentCum = 0e0;
	  momentSqrdCum = 0e0;
	  
	  for (int run = 0; run < L; run++)
	    {
	      momentCum += momentsAllRunsInterp[run][mom][timeStep];
	      momentSqrdCum += momentsAllRunsInterp[run][mom][timeStep]*momentsAllRunsInterp[run][mom][timeStep];
	    }
	  momMean = momentCum/L;
	  momVar = momentSqrdCum/L - momMean*momMean;
	  momentsFile << momMean << "\t" << momVar << "\t";
	}
      momentsFile << endl;  
      
    }
  cout << "Done!" << endl;
  
  // Let's save the state of the random number generator
  if(saveRandState)
    {
      saveState(mtrand);
    }
  
  // Close files
  diamsFile.close();
  momentsFile.close();   
  outputFile.close();
  
  cout << "Simulation complete!" << endl << endl;
  
  return 0;
  
}



  
  
  /*
  cout << "Dumping interpolated moments"<< endl;
  
  // Let's dump the interpolated moments from the first run - works
  for (int run = 0; run < L; run++)
    {
      // output initial condition
      for(int elem = 0; elem < noMoments+1; elem++)
	cout << momentsAllRuns[run][elem][0]<< "\t";
      cout << endl;
      
      
      for (int step = 0; step < momentsAllRunsInterp[run][0].size();step++) // count the # time steps in this run
	{
	  for(int elem = 0; elem < noMoments+1; elem++)
	    {
	      cout << momentsAllRunsInterp[run][elem][step] << "\t";
	    }
	  cout << endl;
	}
    }
  
  */
  



  

  /*
for(vector<vector<vector<double> > >::iterator iterRun = momentsAllRuns.begin(); iterRun != momentsAllRuns.end(); ++iterRun)
    {
      cout << "iterRun = "  << endl;
      for(vector<double>::iterator iterStep=(*iterRun)[0].begin(); iterStep != (*iterRun)[0].end(); ++iterStep)
	{
	  cout << "timeStep  = " << *timeStep;

	}
    }
  */

  
  /*
  // Let's dump the moments from the first run - works
  for (int run = 0; run < 1; run++)
    {
      for (int step = 0; step < momentsAllRuns[run][0].size();step++) // count the # time steps in this run
	{
	  for(int elem = 0; elem < noMoments+1; elem++)
	    {
	      momentsFile << momentsAllRuns[run][elem][step] << "\t";
	    }
	  momentsFile << endl;
	}
    }
  */
      

  

  /*
  vector<vector<vector<double> > >::iterator iterRun = momentsAllRuns.begin();
  cout << "noElems = " << iterRun->size() << endl;
  cout << "steps in first elem = " << (*iterRun)[0].size() << endl;

  for(int elem = 0; elem < iterRun->size(); elem++)
    {
      for (int step = 0; step < (*iterRun)[elem].size(); step++)
	{
	  cout << (*iterRun)[elem][step] << "\t" <<endl;
	}
      cout << endl;
      
      //
      
      
    }

  */


  /*
  // Let's dump the moments from the first run
  vector<vector<vector<double> > >::iterator iterRun = momentsAllRuns.begin();
  for(vector<vector<double> >::iterator iterElem = iterRun->begin(); iterElem != iterRun->end(); ++iterElem)
    {
      for(vector<double>::iterator iterLine = iterElem->begin(); iterLine != iterElem->end(); ++iterLine)
	{
	  momentsFile << *iterLine << "\t";
	}
      momentsFile << endl;
    }
  */

  /*
  // Dump all info - to screen
  for(vector<vector<vector<double> > >::iterator iterRun = momentsAllRuns.begin(); iterRun != momentsAllRuns.end(); ++iterRun)
    {
      cout << "Run: " << iterRun - momentsAllRuns.begin() + 1 << endl;
      for(vector<vector<double> >::iterator iterElem = iterRun->begin(); iterElem != iterRun->end(); ++iterElem)
	{
	  cout << "Element: " << iterElem -  iterRun->begin() << "\t";
	  for(vector<double>::iterator iterLine = iterElem->begin(); iterLine != iterElem->end(); ++iterLine)
	    {
	      cout << "Timestep = " << iterLine -  iterElem->begin() << "value = " << *iterLine << "\t";
	    }
	  cout << endl;
	}

    }

  */



