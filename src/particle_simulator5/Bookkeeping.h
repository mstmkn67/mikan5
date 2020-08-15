#ifndef _BOOKKKEPING_H_
#define _BOOKKEEPING_H_

#include "Particle.h"
#include "BoundaryCondition.h"
#include <map>
#include <vector>
using namespace std;

class Bookkeeping{
public:
	Bookkeeping(map<int,vector<Particle*> >* particle);
	virtual ~Bookkeeping();
	virtual void update();
	//virtual void interaction_table_clear();
	vector<Particle*> particle;//particles
	map<Bead*,vector<Bead*> > interaction_table;
	BoundaryCondition* bc;
protected:
	virtual void bookParticle(Particle* p0,Particle* p1);
private:
};

class BeadsBookkeeping:public Bookkeeping{
public:
	BeadsBookkeeping(map<int,vector<Particle*> >* particle,double radius);
	virtual ~BeadsBookkeeping();
	//virtual void update();
private:
	virtual void bookParticle(Particle* p0,Particle* p1);
	double radius2;
};

class ParticlesBookkeeping:public Bookkeeping{
public:
	ParticlesBookkeeping(map<int,vector<Particle*> >* particle,double radius);
	virtual ~ParticlesBookkeeping();
	//virtual void update();
protected:
	virtual void bookParticle(Particle* p0,Particle* p1);
private:
	double radius2;
};

#endif // _BOOKKEEPING_H_
