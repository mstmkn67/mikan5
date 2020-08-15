//Force

#ifndef _FORCE_H_
#define _FORCE_H_

#include "Particle.h"

class Force{
public:
	Force(const vector<Particle*>& p):particles(p){};
	virtual ~Force(){};
	virtual void update()=0;
protected:
	vector<Particle*> particles;
};

#endif // _FORCE_H_
