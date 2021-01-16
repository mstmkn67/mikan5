#ifndef RNG_H
#define RNG_H
//   ___________
// _| class RNG |_____________________________________________________________
// This is the Abstract Base Class for all Random Number Generators
// It contains the absolute minimum that any RNG should contain:
// 	Seed	- seeds the RNG
//	lRand	- returns a uint32 (unsigned long)
//  dRand	- returns a double in the range [0.0,1.0)
// It also contains an ID function:
//	Info	- returns info about the RNG (author, version, etc)
//

// first we define the info returned by the Info member
struct RNGInfoItem
{
	const char *parameter;	// "Author", "Version", etc
	const char *value;
};

#include <stdio.h>

#ifdef __BEOS__
#include <SupportDefs.h>
#else
// non-BeOs systems need to ensure that this is a 32-bit int!
typedef unsigned long uint32;
#endif

struct RNGInfo
{
	unsigned nitems;
	RNGInfoItem *item;

	void PrintToStream (FILE *fp = stderr) const;
} ;

// and now the Base Class
//	- the pure virtuals make this an abstract base class

class RNG
{
	static RNGInfo *info;
public:
	virtual void Seed (uint32) = 0;
	virtual uint32 lRand() = 0;
	virtual double dRand();

	virtual const RNGInfo *Info() { return info; };
} ;

#endif
