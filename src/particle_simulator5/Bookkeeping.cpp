#include "BookKeeping.h"

//
Bookkeeping::Bookkeeping(map<int,vector<Particle*> >* p){
	particle.clear();
	map<int,vector<Particle*> >::iterator i=p->begin();
	for(;i!=p->end();i++){
		//vector<Particle*>::iterator j=i->second.be
		particle.insert(particle.end(),i->second.begin(),i->second.end());
	}
}

Bookkeeping::~Bookkeeping(){}

void Bookkeeping::update(){
	interaction_table.clear();
	int a=0;
	vector<Particle*>::iterator i=particle.begin();
	for(;i!=particle.end();i++,a++){
		int b=0;
		vector<Particle*>::iterator j=particle.begin();
		for(;j!=i;j++,b++){
			bookParticle(*i,*j);
		}
	}
}

void Bookkeeping::bookParticle(Particle* i,Particle* j){
	vector<Bead*>::iterator a=i->bead.begin();
	for(;a!=i->bead.end();a++){
		vector<Bead*>::iterator b=j->bead.begin();
		for(;b!=j->bead.end();b++){
			interaction_table[*a].push_back(*b);
		}
	}
}

//
BeadsBookkeeping::BeadsBookkeeping(map<int,vector<Particle*> >* p,double r)
:radius2(r*r),Bookkeeping(p){}

BeadsBookkeeping::~BeadsBookkeeping(){}

void BeadsBookkeeping::bookParticle(Particle* i,Particle* j){
	Vector3d r=i->position-j->position;
	if(r.length2()>radius2)return;
	vector<Bead*>::iterator a=i->bead.begin();
	for(;a!=i->bead.end();a++){
		vector<Bead*>::iterator b=j->bead.begin();
		for(;b!=j->bead.end();b++){
			interaction_table[*a].push_back(*b);
		}
	}
}

//
ParticlesBookkeeping::ParticlesBookkeeping(map<int,vector<Particle*> >* p,double r)
:radius2(r*r),Bookkeeping(p){}

ParticlesBookkeeping::~ParticlesBookkeeping(){}

void ParticlesBookkeeping::bookParticle(Particle* i,Particle* j){
	vector<Bead*>::iterator a=i->bead.begin();
	for(;a!=i->bead.end();a++){
		vector<Bead*>::iterator b=j->bead.begin();
		for(;b!=j->bead.end();b++){
			Vector3d r=(*a)->position-(*b)->position;
			if(r.length2()>radius2)continue;
			interaction_table[*a].push_back(*b);
		}
	}
}

