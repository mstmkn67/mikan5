#ifndef MERSENNE_TWISTER
#define MERSENNE_TWISTER

// Header file for MT19937:
// lRand() generates one pseudorandom unsigned integer (32bit)
//  which is uniformly distributed among 0 to 2^32-1
// dRand() generates one pseudorandom real number
//  which is uniformly distributed among 0.0 to <1.0
// Coded by Takuji Nishimura, considering the suggestions by
// Topher Cooper and Marc Rieffel in July-Aug. 1997.
// Modified by John Cooke 1998-04-20
//  to use standard Seed(), lRand(), dRand() functions

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Library General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later
// version.
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU Library General Public License for more details.
// You should have received a copy of the GNU Library General
// Public License along with this library; if not, write to the
// Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
// 02111-1307  USA

// Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
// When you use this, send an email to: matumoto@math.keio.ac.jp
// with an appropriate reference to your work.

#ifndef RNG_H
#include "rng.h"
#endif

// Period parameters
#define MT_BUFSZ 624

class Twister : public RNG
{
private:
    uint32   state[MT_BUFSZ+1];  // state vector + 1 extra to not violate ANSI C
    uint32   *next;          // next random value is computed from here
    int      left;      // can *next++ this many times before reloading

    static RNGInfo *info;

    uint32 reloadMT(void);

public:
// we must override these two pure virtual functions ...
    void Seed (uint32);
    uint32 lRand ();

// we override these two ...
    double dRand();
    const RNGInfo *Info() { return info; };

// we will have a constructor too
    Twister (uint32);
} ;

#endif
