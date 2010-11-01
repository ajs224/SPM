// --------------------------------------
// Random numbers
// A. J. Smith (ajs224@cam.ac.uk)
//
//---------------------------------------


#include <iostream>
#include <fstream>
#include "random.h"

using namespace std;

namespace ajsRandom
{
  int getSeed()
  {
    ifstream randSeed("/dev/urandom");
    char tmp[sizeof(int)];
    randSeed.read(tmp,sizeof(int));
    randSeed.close();
    int * seed = reinterpret_cast<int *>(tmp);
    return (*seed);
  }

  // Declare a Mersenne Twister random number generator
  // Either reload the seeds or generate random new ones
  bool loadRandState=true;
  bool saveRandState=false;

  MTRand myRand(bool loadRandState)
  {
    MTRand mtrand;
    if(loadRandState)
      {
	cout << "Loading previous state of random number generator..." <<endl;
	
	// Read state from file
	ifstream stateIn( "state.dat" );
	if( stateIn )
	  {
	    stateIn >> mtrand;
	    stateIn.close();
	  }
      }
    else
      {
	cout << "Generating new seeds for random number generator..." <<endl;
	// We seed with an array of values rather than a single integer
	// to access the full 2^19937-1 possible sequences.
	MTRand::uint32 seed[ MTRand::N ];
	for( int n = 0; n < MTRand::N; ++n )
	  {
	    //seed[n] = 23 * n;  // fill with anything
	    seed[n] = getSeed();  // fill with anything
	  }
	
	MTRand mtrand( seed );
      }
    
    //cout << "Done" <<endl;
     
    
    //double m1 = mtrand();
    //cout << "rand no =" << m1 << endl;

    return mtrand;

  }

  void saveState(MTRand mtrand)
  {
    cout << "Saving state of Random number generator..." << endl;
    
    // Save stream to a file
    ofstream stateOut( "state.dat" );
    if( stateOut )
      {
	stateOut << mtrand;
	stateOut.close();
      }
  }
  
}//namespace  ajsRandom
