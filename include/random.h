#ifndef RANDOM_H
#define RANDOM_H

#include "MersenneTwister.h"

namespace ajsRandom
{

  extern int getSeed();

  extern MTRand myRand(bool); // Mersenne twister
  extern void saveState(MTRand);
  extern bool loadRandState; // load old random state?
  extern bool saveRandState; // save current random state?
  
} //namespace  ajsRandom

#endif
