// RNG info
#include "rng.h"

//   
// _||
// This is the static info for the RNG class
//
static RNGInfoItem items0[] =
{
	{ "Class", "RNG"},
	{ "Description", "Abstract Base Class for Random Number Generators" },
	{ "Version", "1.0 : 1998-04-19" },
	{ "Author", "John Cooke" },
	{ "Email", "bebox@dircon.co.uk" }
} ;

static RNGInfo info0 = { 5, items0 } ;

RNGInfo *RNG::info = &info0;

//   
//_| RNG MEMBER FUNCTIONS |

// dRand calls lRand and returns a double within the range 0.0-1.0
// to achieve this we can divide by 2**32
// or, which should be faster, multiply by 1/(2**32)
//
#define L2D (1.0/4294967296.0)

double RNG::dRand()
{
	double d = (double)lRand();
	return d * L2D;
}

//   
// _| RNGInfo MEMBER FUNCTIONS |
//
#include <stdio.h>

void RNGInfo::PrintToStream (FILE *fp) const
{
	const RNGInfoItem *it = item;

	for(unsigned i=0; i<nitems; i++,it++)
		fprintf(fp, "%s: %s\n", it->parameter, it->value);
}

