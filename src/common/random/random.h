#ifndef _RANDOM_H_INCLUDED_
#define _RANDOM_H_INCLUDED_
// $Id: random.h,v 1.2 2003/05/23 16:28:38 takimoto Exp $
/*====================================================================

  random.h - Easy interface for random number generators

  UniformRandomNumber:  Generate uniform random numbers  
  GaussianRandomNumber: Generate gaussian random numbers
  
  ver 1.0	1998	Masao Doi
  Define interface.
  Use the algorithm in "Numerical Recepie in C".

  version 2.0 1999/4/14 T.Aoyagi
  Use "Mersenne Twister" algorithm.
  rng.h, rng.cpp, twister.h, twister.cpp -- by John Cooke.
  
  version 2.1 1999/5/21 T.Aoyagi
  Get seed from time().

  version 2.2 2001/9/26 J.Takimoto
  Some clean up.

	Usage:
	[1] Construct.

			UniformRandomNumber	uran;

		seed is generated from current time.
	
	[2] Set parameters, if necessary.

			uran.setMax(-10.0).setMin(-10.0);

		default is max=1, min=0.

	[3] Use operator() to get random numbers.

			x = uran();

 ====================================================================*/

#include "rng.h"

//----- UniformRandomNumber -----

class UniformRandomNumber {
public:
	explicit UniformRandomNumber(uint32 seed=0);
	~UniformRandomNumber() { delete rng; }
	double operator() ();
	UniformRandomNumber&  setMax(double max_=1.0);
	UniformRandomNumber&  setMin(double min_=0.0);
	UniformRandomNumber&  setSeed(uint32 seed);
private:
	RNG* rng;
	double max;
	double min;
	double delta;	// = max - min
};

inline double
UniformRandomNumber::operator() ()
{
	return rng->dRand()*delta + min;
}

//----- GaussianRandomNumber -----

class GaussianRandomNumber {
public:
	explicit GaussianRandomNumber(uint32 seed=0);
	~GaussianRandomNumber(){delete rng;}
	double operator() ();
	GaussianRandomNumber&  setSigma(double sigma_=1.0);
	GaussianRandomNumber&  setMean(double mean_=0.0);
	GaussianRandomNumber&  setSeed(uint32 seed);
private:
	RNG* rng;
	double sigma;
	double mean;
	bool use_saved;
	double saved;
};

#endif // _RANDOM_H_INCLUDED_
