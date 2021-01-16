// $Id: random.cpp,v 1.2 2003/05/23 16:28:38 takimoto Exp $
//----------------
//	random.cpp
//----------------
#include <math.h>
#include <time.h>
#include "random.h"
#include "twister.h"

static uint32 default_seed = time((time_t*)0);

// if seed==0, use defualt_seed

static uint32
myseed(uint32 seed)
{
	if( seed == 0 ) {
		seed = default_seed;
		default_seed += 2;	// why 2? see twister.cpp
	}
	return seed;
}

//----- UniformRandomNumber -----

//UniformRandomNumber::UniformRandomNumber(uint32 seed=0)
UniformRandomNumber::UniformRandomNumber(uint32 seed)
{
	rng = new Twister(myseed(seed));

	max = 1.0;
	min = 0.0;
	delta = 1.0;
}

UniformRandomNumber&
UniformRandomNumber::setMax(double max_) {
	max = max_; 
	delta = max - min;
	return *this;
}

UniformRandomNumber&
UniformRandomNumber::setMin(double min_){
	min = min_; 
	delta = max - min;
	return *this;
}

UniformRandomNumber&
UniformRandomNumber::setSeed(uint32 seed){
	rng->Seed(seed);
	return *this;
}

//----- GaussianRandomNumber -----

//GaussianRandomNumber::GaussianRandomNumber(uint32 seed=0)
GaussianRandomNumber::GaussianRandomNumber(uint32 seed)
{  
	rng = new Twister(myseed(seed));
	
	sigma = 1.0;
	mean = 0.0;
	use_saved = false;
}

GaussianRandomNumber&  
GaussianRandomNumber::setSigma(double sigma_)
{
	sigma = sigma_;
	return *this;
}

GaussianRandomNumber&  
GaussianRandomNumber::setMean(double mean_)
{
	mean = mean_;
	return *this;
}

GaussianRandomNumber&  
GaussianRandomNumber::setSeed(uint32 seed)
{
	rng->Seed(seed);
	return *this;
}

//--- Box-Muller method

double
GaussianRandomNumber::operator() ()
{
	double x, y, rr, g;

	if( use_saved ) {
		use_saved = false;
		return saved;
	}
	else {
		do {		// generate a point within a unit circle
			x = 2.0*rng->dRand() - 1.0;
			y = 2.0*rng->dRand() - 1.0;
			rr = x*x + y*y;
		} while (rr >= 1.0 || rr == 0.0);

		g = sigma*sqrt(-2.0*log(rr)/rr);
		saved = y*g + mean;
		use_saved = true;
		return x*g + mean ;
	}
}
